#!/usr/bin/env python3
"""Emit a machine-readable dataset record for a processed study.

The OriDB data hub needs to describe each ingested dataset in a uniform way: what was
measured, in which strain, from which archived runs, against which assembly, and by
which pipeline. This collects that from the samplesheet, the comparison sheet and the
pipeline's own output statistics, so the record cannot drift from the data it describes.
"""
import argparse, csv, datetime, hashlib, json, os, subprocess


def read_kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                out[k] = v
    return out


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def git_revision(repo):
    try:
        return subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--results", required=True, help="pipeline --outdir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dataset-id", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--doi", default="")
    ap.add_argument("--pubmed", default="")
    ap.add_argument("--series", default="")
    ap.add_argument("--assembly", default="sacCer3")
    ap.add_argument("--organism", default="Saccharomyces cerevisiae")
    ap.add_argument("--bin-size", type=int, default=1000)
    ap.add_argument("--repo", default=".")
    a = ap.parse_args()

    samples = list(csv.DictReader(open(a.samplesheet)))
    comps = list(csv.DictReader(open(a.comparisons)))

    sample_records = []
    for s in samples:
        counts = os.path.join(a.results, "counts",
                              f"{s['sample']}.{a.bin_size}bp.stats.txt")
        qc = os.path.join(a.results, "qc", f"{s['sample']}.qc.txt")
        rec = {k: s[k] for k in
               ("sample", "run", "platform", "layout", "variant", "role", "strain",
                "timepoint", "series") if k in s}
        rec["counts"] = read_kv(counts)
        rec["qc"] = read_kv(qc)
        sample_records.append(rec)

    track_records = []
    for c in comps:
        stats = read_kv(os.path.join(a.results, "profiles",
                                     f"{c['comparison']}.normstats.txt"))
        rec = {
            "comparison": c["comparison"],
            "replicating_sample": c["rep"],
            "non_replicating_sample": c["nonRep"],
            "variant": c["variant"],
            "normalisation_mode": c["norm_mode"],
            "description": c["description"],
            "statistics": stats,
            "files": {},
        }
        for kind, rel in (("ratio_bed", f"profiles/{c['comparison']}.ratio.bed"),
                          ("smoothed_bed", f"profiles/{c['comparison']}.smooth.bed"),
                          ("ratio_table", f"profiles/{c['comparison']}.ratio.tsv"),
                          ("bigwig", f"bigwig/{c['comparison']}.smooth.bw")):
            p = os.path.join(a.results, rel)
            if os.path.exists(p):
                rec["files"][kind] = {"path": rel,
                                      "bytes": os.path.getsize(p),
                                      "md5": md5(p)}
        track_records.append(rec)

    record = {
        "dataset_id": a.dataset_id,
        "title": a.title,
        "organism": a.organism,
        "assembly": a.assembly,
        "bin_size_bp": a.bin_size,
        "measurement": "relative DNA copy number (replication timing)",
        "scale": "1 = unreplicated, 2 = fully replicated",
        "source": {
            "doi": a.doi,
            "pubmed": a.pubmed,
            "archive_series": [s for s in a.series.split(",") if s],
            "runs": sorted({s["run"] for s in samples if s.get("run")}),
        },
        "provenance": {
            "pipeline": "Muller-replication-timing",
            "pipeline_revision": git_revision(a.repo),
            "generated": datetime.datetime.now(datetime.timezone.utc)
                                  .strftime("%Y-%m-%dT%H:%M:%SZ"),
        },
        "samples": sample_records,
        "tracks": track_records,
    }

    with open(a.out, "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    print(f"{a.out}: {len(sample_records)} samples, {len(track_records)} tracks")


if __name__ == "__main__":
    main()
