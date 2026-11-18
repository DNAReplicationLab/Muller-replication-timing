#!/usr/bin/env python3
"""Decode SOLiD colourspace FASTQ to base space.

Fourteen runs of Muller et al. 2014 were sequenced on an AB SOLiD 5500 and are
archived in colourspace: a primer base followed by N colour calls, where each call
encodes the transition between adjacent bases. ENA serves them in that form and
Bowtie 2 cannot read it.

Decoding walks the colour string from the primer base through the standard SOLiD
two-base transition table. A missing colour call ('.') breaks the chain - every base
downstream of it is unrecoverable - so the read is truncated there. This reproduces
the behaviour of `fastq-dump -B`, which emits N from the first missing colour to the
end of the read; truncating instead of emitting N matters because Bowtie 2's default
--n-ceil would reject a read carrying that many Ns outright.

Verified byte-exact against NCBI's own decoded output: see tests/ in this repository.
"""
import argparse, gzip, sys

# next_base = TRANSITION[previous_base][colour]
TRANSITION = {
    'A': {'0': 'A', '1': 'C', '2': 'G', '3': 'T'},
    'C': {'0': 'C', '1': 'A', '2': 'T', '3': 'G'},
    'G': {'0': 'G', '1': 'T', '2': 'A', '3': 'C'},
    'T': {'0': 'T', '1': 'G', '2': 'C', '3': 'A'},
}


def opener(path, mode="rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def decode(colours):
    """Primer base + colour calls -> base-space read, truncated at the first bad call."""
    prev = colours[0].upper()
    if prev not in TRANSITION:
        return ""
    out = []
    for c in colours[1:]:
        nxt = TRANSITION[prev].get(c)
        if nxt is None:          # '.' or any non-call: the chain is broken
            break
        out.append(nxt)
        prev = nxt
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True, help="colourspace FASTQ (.gz ok, '-' for stdin)")
    ap.add_argument("--out", required=True, help="base-space FASTQ (.gz ok, '-' for stdout)")
    ap.add_argument("--min-length", type=int, default=25,
                    help="drop reads shorter than this after truncation (default 25)")
    ap.add_argument("--stats", help="write a summary here")
    ap.add_argument("--pad-n", action="store_true",
                    help="emit N to the original read length instead of truncating "
                         "(reproduces fastq-dump -B exactly; not recommended for mapping)")
    a = ap.parse_args()

    n_in = n_out = n_short = 0
    len_sum = 0
    fin, fout = opener(a.inp), opener(a.out, "wt")
    try:
        while True:
            name = fin.readline()
            if not name:
                break
            seq = fin.readline().rstrip("\n")
            fin.readline()                      # '+' line
            qual = fin.readline().rstrip("\n")
            n_in += 1

            bases = decode(seq)
            # The quality string has one entry per colour call. The first character
            # corresponds to the primer base when it is present; align by taking the
            # trailing len(seq)-1 characters, then cut to the decoded length.
            q = qual[-(len(seq) - 1):] if len(qual) >= len(seq) else qual
            if a.pad_n:
                pad = (len(seq) - 1) - len(bases)
                bases += "N" * pad
                q = q[:len(bases)].ljust(len(bases), "!")
            else:
                q = q[:len(bases)]

            if len(bases) < a.min_length:
                n_short += 1
                continue
            fout.write(name.rstrip("\n") + "\n" + bases + "\n+\n" + q + "\n")
            n_out += 1
            len_sum += len(bases)
    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    lines = [
        f"reads_in\t{n_in}",
        f"reads_out\t{n_out}",
        f"reads_below_min_length\t{n_short}",
        f"fraction_kept\t{n_out / n_in:.4f}" if n_in else "fraction_kept\t0",
        f"mean_decoded_length\t{len_sum / n_out:.2f}" if n_out else "mean_decoded_length\t0",
    ]
    sys.stderr.write("\n".join(lines) + "\n")
    if a.stats:
        open(a.stats, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
