import argparse
import logging

from . import __version__
from .quant import quantify


def read_lengths(spec):
    try:
        lo, hi = (int(x) for x in spec.split("-"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected LO-HI, got {spec!r}")
    if lo > hi:
        raise argparse.ArgumentTypeError(f"LO > HI in {spec!r}")
    return lo, hi


def main(argv=None):
    p = argparse.ArgumentParser(prog="ribokit")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("quant", help="P-site offsets and CDS quantification from a transcriptome BAM")
    q.add_argument("--bam", required=True, help="transcriptome alignments (reference names = transcript ids)")
    q.add_argument("--gtf", required=True)
    q.add_argument("--fasta", required=True, help="genome FASTA of the GTF")
    q.add_argument("--read-lengths", required=True, type=read_lengths, metavar="LO-HI")
    q.add_argument("--out-prefix", required=True,
                   help="writes <prefix>.quant.tsv, .offsets.tsv, .ties.tsv, .stats.tsv, .psites.tsv")
    q.add_argument("--offsets", help="use this offsets.tsv instead of estimating offsets from the BAM")
    q.add_argument("--min-offset-support", type=float, default=30,
                   help="reads spanning a start or stop codon needed to give a read length offsets")
    q.add_argument("--tol", type=float, default=1e-3,
                   help="EM stops when one step moves no expected count by this many reads")
    q.add_argument("--max-iter", type=int, default=100_000,
                   help="give up after this many EM steps (SQUAREM counts 2-3 per cycle)")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    quantify(a.bam, a.gtf, a.fasta, a.read_lengths, a.out_prefix, offsets_path=a.offsets,
             min_support=a.min_offset_support, tol=a.tol, max_iter=a.max_iter)


if __name__ == "__main__":
    main()
