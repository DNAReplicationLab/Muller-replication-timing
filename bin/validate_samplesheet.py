#!/usr/bin/env python3
"""Check a samplesheet and comparison sheet are internally consistent.

Run by CI, and worth running by hand before a long execution: every failure it
catches would otherwise appear hours into a pipeline run.
"""
import argparse, csv, re, sys

SAMPLE_COLS = {"sample", "run", "platform", "layout", "variant", "role",
               "strain", "timepoint", "series"}
COMP_COLS = {"comparison", "rep", "nonRep", "variant", "norm_mode",
             "norm_value", "description"}
PLATFORMS = {"ILLUMINA", "ABI_SOLID"}
LAYOUTS = {"single", "paired"}
ROLES = {"rep", "nonRep"}
VARIANTS = {"sortseq", "syncseq", "mfaseq", "hu"}
NORM_MODES = {"auto", "rFactor", "upperLimit"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--comparisons", required=True)
    a = ap.parse_args()
    errors = []

    def ragged(path, n_expected):
        """Rows with more fields than the header - usually an unquoted comma in a
        description. Nextflow's splitCsv aborts the whole run on these with a message
        that does not say which file is at fault."""
        bad = []
        with open(path) as fh:
            for i, row in enumerate(csv.reader(fh), 1):
                if row and len(row) != n_expected:
                    bad.append((i, len(row)))
        return bad

    samples = list(csv.DictReader(open(a.samplesheet)))
    if not samples:
        sys.exit(f"{a.samplesheet} is empty")
    missing = SAMPLE_COLS - set(samples[0])
    if missing:
        errors.append(f"samplesheet missing columns: {sorted(missing)}")

    for path, n in ((a.samplesheet, len(SAMPLE_COLS)), (a.comparisons, len(COMP_COLS))):
        for line_no, width in ragged(path, n):
            errors.append(f"{path} line {line_no}: {width} fields, header has {n} "
                          f"- an unquoted comma in a field?")

    seen = set()
    for i, r in enumerate(samples, 2):
        where = f"samplesheet line {i} ({r.get('sample')})"
        if r["sample"] in seen:
            errors.append(f"{where}: duplicate sample name")
        seen.add(r["sample"])
        if r["platform"] not in PLATFORMS:
            errors.append(f"{where}: platform {r['platform']!r} not in {sorted(PLATFORMS)}")
        if r["layout"] not in LAYOUTS:
            errors.append(f"{where}: layout {r['layout']!r} not in {sorted(LAYOUTS)}")
        if r["role"] not in ROLES:
            errors.append(f"{where}: role {r['role']!r} not in {sorted(ROLES)}")
        if r["variant"] not in VARIANTS:
            errors.append(f"{where}: variant {r['variant']!r} not in {sorted(VARIANTS)}")
        if r["run"] and not re.fullmatch(r"[SED]RR\d+", r["run"]):
            errors.append(f"{where}: {r['run']!r} is not an SRA run accession")

    comps = list(csv.DictReader(open(a.comparisons)))
    missing = COMP_COLS - set(comps[0] if comps else {})
    if missing:
        errors.append(f"comparisons missing columns: {sorted(missing)}")

    by_name = {r["sample"]: r for r in samples}
    for i, c in enumerate(comps, 2):
        where = f"comparisons line {i} ({c.get('comparison')})"
        for role_col, expect in (("rep", "rep"), ("nonRep", "nonRep")):
            name = c[role_col]
            if name not in by_name:
                errors.append(f"{where}: {role_col} {name!r} is not in the samplesheet")
            elif by_name[name]["role"] != expect:
                errors.append(f"{where}: {name!r} has role "
                              f"{by_name[name]['role']!r}, expected {expect!r}")
        if c["norm_mode"] not in NORM_MODES:
            errors.append(f"{where}: norm_mode {c['norm_mode']!r} not in {sorted(NORM_MODES)}")
        if c["norm_mode"] in {"rFactor", "upperLimit"} and not c["norm_value"]:
            errors.append(f"{where}: norm_mode {c['norm_mode']} requires a norm_value")
        if c["norm_mode"] == "rFactor" and c["norm_value"]:
            v = float(c["norm_value"])
            if not 1.0 < v <= 2.0:
                errors.append(f"{where}: rFactor {v} outside the 1-2 copy number scale")
        if c["variant"] == "syncseq" and not re.search(r"\d+min", c["comparison"]):
            errors.append(f"{where}: sync-seq comparison names must embed the "
                          f"timepoint as e.g. '35min' - Trep ordering depends on it")

    if errors:
        print("\n".join(f"ERROR  {e}" for e in errors))
        sys.exit(1)
    print(f"OK: {len(samples)} samples, {len(comps)} comparisons, no inconsistencies")


if __name__ == "__main__":
    main()
