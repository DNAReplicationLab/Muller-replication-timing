#!/usr/bin/env python3
"""Convert a 5-column BED of binned values to bigWig.

UCSC's bedGraphToBigWig is not reliably available on all platforms (the bioconda
osx-arm64 build links against a Homebrew libpng), so this uses pyBigWig, which writes
a correct file including its total-summary block.
"""
import argparse, sys
import pyBigWig


BIGWIG_MAGIC = 0x888FFC26
TOTAL_SUMMARY_OFFSET_POS = 44        # bbiFile header: where totalSummaryOffset lives


def bigwig_summary(path):
    """Read a bigWig's total-summary block directly.

    pyBigWig 0.3.30's header() returns minVal and maxVal cast to int, so a profile
    spanning 0.933 to 2.276 is reported as 0 to 2. The file itself is correct - only
    that accessor is lossy - so anything that checks a track's value range has to read
    the summary block rather than trust header().
    """
    import struct
    with open(path, "rb") as fh:
        head = fh.read(64)
        if struct.unpack("<I", head[:4])[0] == BIGWIG_MAGIC:
            endian = "<"
        elif struct.unpack(">I", head[:4])[0] == BIGWIG_MAGIC:
            endian = ">"
        else:
            raise ValueError(f"{path} is not a bigWig file")
        offset = struct.unpack(endian + "Q",
                               head[TOTAL_SUMMARY_OFFSET_POS:TOTAL_SUMMARY_OFFSET_POS + 8])[0]
        if offset == 0:
            raise ValueError(f"{path} has no total-summary block")
        fh.seek(offset)
        valid, vmin, vmax, total, total_sq = struct.unpack(endian + "Qdddd", fh.read(40))
    return {"validCount": valid, "minVal": vmin, "maxVal": vmax,
            "sumData": total, "sumSquared": total_sq,
            "mean": total / valid if valid else float("nan")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sizes = {}
    for line in open(a.chrom_sizes):
        c, n = line.split()[:2]
        sizes[c] = int(n)

    rows = []
    for line in open(a.bed):
        if not line.strip() or line.startswith(("track", "#", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        chrom, start, end, value = f[0], int(f[1]), int(f[2]), float(f[4])
        if chrom not in sizes:
            continue
        end = min(end, sizes[chrom])
        if end <= start:
            continue
        rows.append((chrom, start, end, value))

    if not rows:
        sys.exit("no usable intervals in %s" % a.bed)

    order = list(sizes)
    rows.sort(key=lambda r: (order.index(r[0]), r[1]))

    bw = pyBigWig.open(a.out, "w")
    bw.addHeader([(c, sizes[c]) for c in order if any(r[0] == c for r in rows)])
    bw.addEntries([r[0] for r in rows], [r[1] for r in rows],
                  ends=[r[2] for r in rows], values=[r[3] for r in rows])
    bw.close()

    vmin = min(r[3] for r in rows)
    vmax = max(r[3] for r in rows)
    summary = bigwig_summary(a.out)
    if abs(summary["minVal"] - vmin) > 1e-5 or abs(summary["maxVal"] - vmax) > 1e-5:
        sys.exit(f"{a.out}: written summary is {summary['minVal']}-{summary['maxVal']}, "
                 f"data are {vmin}-{vmax}")

    print(f"{a.out}: {len(rows)} intervals, {len({r[0] for r in rows})} sequences, "
          f"value range {vmin:.4f}-{vmax:.4f} (summary block verified)")


if __name__ == "__main__":
    main()
