#!/usr/bin/env python3
"""Correlate two binned profiles on their shared bins.

Used to check a newly generated profile against a published one, and to compare
replicates. Bins are matched on (chrom, start, end), so both files must use the same
assembly and bin size.
"""
import argparse, csv, sys
from statistics import mean


def load(path, value_col=None):
    """Read a BED-like file of binned values.

    The value is taken from the last column unless --value-col says otherwise, so
    both BED4 (chrom,start,end,value) and BED5 (with a name column) work without
    the caller having to know which.
    """
    out = {}
    for line in open(path):
        if not line.strip() or line.startswith(("track", "#", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        col = value_col if value_col is not None else len(f) - 1
        try:
            out[(f[0], int(f[1]), int(f[2]))] = float(f[col])
        except (ValueError, IndexError):
            continue
    if not out:
        sys.exit(f"no usable rows read from {path}")
    return out


def pearson(x, y):
    n = len(x)
    mx, my = mean(x), mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / (sxx * syy) ** 0.5 if sxx and syy else float("nan")


def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):          # average ranks within ties
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="BED5 profile A")
    ap.add_argument("--b", required=True, help="BED5 profile B")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--exclude-chrom", default="chrM")
    ap.add_argument("--value-col", type=int, default=None,
                    help="0-based column holding the value (default: last column)")
    ap.add_argument("--out")
    args = ap.parse_args()

    A, B = load(args.a, args.value_col), load(args.b, args.value_col)
    excl = set(args.exclude_chrom.split(",")) if args.exclude_chrom else set()
    keys = sorted(k for k in set(A) & set(B) if k[0] not in excl)
    if not keys:
        sys.exit("no shared bins - different assembly or bin size?")

    x = [A[k] for k in keys]
    y = [B[k] for k in keys]
    rx, ry = rank(x), rank(y)

    # slope of B on A through the origin-free least squares fit, which is what
    # recovers a scale factor between two differently normalised profiles
    mx, my = mean(x), mean(y)
    sxx = sum((a - mx) ** 2 for a in x)
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx if sxx else float("nan")

    res = {
        "profile_a": args.label_a,
        "profile_b": args.label_b,
        "bins_a": len(A),
        "bins_b": len(B),
        "bins_shared": len(keys),
        "pearson_r": round(pearson(x, y), 5),
        "spearman_rho": round(pearson(rx, ry), 5),
        "mean_a": round(mx, 5),
        "mean_b": round(my, 5),
        "slope_b_on_a": round(slope, 5),
        "intercept": round(my - slope * mx, 5),
    }
    for k, v in res.items():
        print(f"{k}\t{v}")
    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(res))
            w.writeheader()
            w.writerow(res)


if __name__ == "__main__":
    main()
