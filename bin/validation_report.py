#!/usr/bin/env python3
"""Correlate a set of new profiles against reference profiles and tabulate the result.

Takes a manifest so the same driver works for any dataset: validating reprocessed
profiles against previously published ones, or comparing replicates against each other.

Manifest is TSV with a header and the columns:
    label    new_bed    reference_bed    [note]
"""
import argparse, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_profiles import load, pearson, rank      # noqa: E402

FIELDS = ["label", "bins_new", "bins_reference", "bins_shared",
          "pearson_r", "spearman_rho", "mean_new", "mean_reference",
          "slope", "intercept", "note"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude-chrom", default="chrM")
    ap.add_argument("--min-r", type=float, default=None,
                    help="exit non-zero if any comparison falls below this correlation")
    a = ap.parse_args()

    excl = set(a.exclude_chrom.split(",")) if a.exclude_chrom else set()
    rows = []
    for m in csv.DictReader(open(a.manifest), delimiter="\t"):
        new, ref = load(m["new_bed"]), load(m["reference_bed"])
        keys = sorted(k for k in set(new) & set(ref) if k[0] not in excl)
        if not keys:
            sys.exit(f"{m['label']}: no shared bins between "
                     f"{m['new_bed']} and {m['reference_bed']}")
        x = [ref[k] for k in keys]
        y = [new[k] for k in keys]
        mx = sum(x) / len(x)
        my = sum(y) / len(y)
        sxx = sum((v - mx) ** 2 for v in x)
        slope = sum((p - mx) * (q - my) for p, q in zip(x, y)) / sxx if sxx else float("nan")
        rows.append({
            "label": m["label"],
            "bins_new": len(new),
            "bins_reference": len(ref),
            "bins_shared": len(keys),
            "pearson_r": round(pearson(x, y), 4),
            "spearman_rho": round(pearson(rank(x), rank(y)), 4),
            "mean_new": round(my, 4),
            "mean_reference": round(mx, 4),
            "slope": round(slope, 4),
            "intercept": round(my - slope * mx, 4),
            "note": m.get("note", ""),
        })
        print(f"{rows[-1]['label']:32s} r={rows[-1]['pearson_r']:.4f}  "
              f"rho={rows[-1]['spearman_rho']:.4f}  n={len(keys)}")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out}: {len(rows)} comparisons")

    if a.min_r is not None:
        bad = [r for r in rows if r["pearson_r"] < a.min_r]
        if bad:
            print(f"\n{len(bad)} comparison(s) below r={a.min_r}: "
                  f"{', '.join(r['label'] for r in bad)}")
            sys.exit(1)


if __name__ == "__main__":
    main()
