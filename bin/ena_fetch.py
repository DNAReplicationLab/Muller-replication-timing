#!/usr/bin/env python3
"""Fetch a sequencing run from ENA as gzipped FASTQ, decoding colourspace if needed.

ENA mirrors every SRA run, serves plain FASTQ over HTTPS and is built for bulk
transfer. File paths are resolved through the ENA portal API rather than constructed,
because the directory layout depends on the accession's checksum suffix.

Transfers use curl with resume and stall detection: these files run to several GB and
a silently stalled connection on a multi-hour run is otherwise indistinguishable from
slow progress. Raw downloads land in a persistent cache directory so an interrupted
run resumes rather than restarts, and every file is checked against the MD5 that ENA
publishes before it is used.
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, urllib.request

PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
          "?accession={run}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,"
          "fastq_bytes,read_count,instrument_platform&format=json")


def resolve(run):
    with urllib.request.urlopen(PORTAL.format(run=run), timeout=120) as fh:
        rows = json.load(fh)
    if not rows:
        sys.exit(f"ENA has no read_run record for {run}")
    row = rows[0]
    urls = [u for u in row["fastq_ftp"].split(";") if u]
    md5s = [m for m in row.get("fastq_md5", "").split(";") if m]
    sizes = [int(b) for b in row.get("fastq_bytes", "").split(";") if b]
    if not urls:
        sys.exit(f"ENA lists no FASTQ files for {run}")
    return ["https://" + u for u in urls], md5s, sizes, row


def md5sum(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest, expect_md5=None, expect_size=None, attempts=4):
    """Download with resume; verify size and MD5; retry from scratch on corruption."""
    # Already here and complete? Verify and skip - curl -C - on a complete file asks
    # for a range past the end and the server answers 416.
    if expect_size and os.path.exists(dest) and os.path.getsize(dest) == expect_size:
        if not expect_md5 or md5sum(dest) == expect_md5:
            sys.stderr.write(f"  cached, complete and verified: {dest}\n")
            return
        sys.stderr.write("  cached file failed its checksum; re-downloading\n")
        os.remove(dest)

    for attempt in range(1, attempts + 1):
        if expect_size and os.path.exists(dest) and os.path.getsize(dest) > expect_size:
            os.remove(dest)                      # overshoot means a bad resume offset
        sys.stderr.write(f"  [{attempt}/{attempts}] {url}\n")
        rc = subprocess.run([
            "curl", "-fsSL", "--continue-at", "-",
            "--retry", "10", "--retry-delay", "5", "--retry-all-errors",
            "--speed-limit", "10240", "--speed-time", "120",   # abort a stalled transfer
            "--connect-timeout", "60",
            "-o", dest, url,
        ]).returncode

        size_ok = (not expect_size) or os.path.getsize(dest) == expect_size
        if rc == 0 and size_ok:
            if not expect_md5:
                return
            got = md5sum(dest)
            if got == expect_md5:
                sys.stderr.write(f"  md5 ok ({got})\n")
                return
            sys.stderr.write(f"  md5 MISMATCH: got {got}, ENA says {expect_md5}\n")
            os.remove(dest)                      # corrupt: start over
        else:
            sys.stderr.write(f"  curl exit {rc}, size "
                             f"{os.path.getsize(dest) if os.path.exists(dest) else 0}"
                             f"/{expect_size}; resuming\n")
    sys.exit(f"failed to download {url} after {attempts} attempts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--paired", action="store_true")
    ap.add_argument("--decode-colourspace", action="store_true")
    ap.add_argument("--min-length", type=int, default=25)
    ap.add_argument("--cache-dir", default=".",
                    help="persistent directory for raw downloads, so a restart resumes")
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--stats")
    a = ap.parse_args()

    urls, md5s, sizes, row = resolve(a.run)
    sys.stderr.write(f"{a.run}: ENA lists {len(urls)} file(s), "
                     f"{row.get('read_count')} reads, platform {row.get('instrument_platform')}\n")

    order = sorted(range(len(urls)), key=lambda i: urls[i])
    if a.paired:
        order = [i for i in order if urls[i].endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        if len(order) != 2:
            sys.exit(f"{a.run}: expected two mate files, ENA lists {len(order)}")
        targets = [f"{a.sample}_1.fastq.gz", f"{a.sample}_2.fastq.gz"]
    else:
        order = [min(order, key=lambda i: len(urls[i]))]
        targets = [f"{a.sample}.fastq.gz"]

    os.makedirs(a.cache_dir, exist_ok=True)
    stats = [f"run\t{a.run}", f"sample\t{a.sample}",
             f"ena_read_count\t{row.get('read_count')}"]

    for i, target in zip(order, targets):
        raw = os.path.join(a.cache_dir, f"raw_{a.run}_{os.path.basename(urls[i])}")
        fetch(urls[i], raw,
              expect_md5=md5s[i] if i < len(md5s) else None,
              expect_size=sizes[i] if i < len(sizes) else None)
        stats.append(f"downloaded_bytes\t{os.path.getsize(raw)}")

        if a.decode_colourspace:
            decoder = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "decode_colourspace.py")
            sys.stderr.write(f"  decoding colourspace -> {target}\n")
            res = subprocess.run(
                [sys.executable, decoder, "--in", raw, "--out", target,
                 "--min-length", str(a.min_length)],
                capture_output=True, text=True)
            sys.stderr.write(res.stderr)
            if res.returncode:
                sys.exit(f"colourspace decoding failed for {a.run}")
            stats += [l for l in res.stderr.strip().split("\n") if "\t" in l]
        else:
            if os.path.exists(target):
                os.remove(target)
            try:
                os.link(raw, target)          # same filesystem: no second copy
            except OSError:
                shutil.copy2(raw, target)

        if not a.keep_raw and a.cache_dir != ".":
            pass        # the cache is the point; leave raw files for resume

    if a.stats:
        open(a.stats, "w").write("\n".join(stats) + "\n")


if __name__ == "__main__":
    main()
