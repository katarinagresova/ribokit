import argparse
import logging

from . import __version__, orfs, score
from .quant import quantify


def read_lengths(spec):
    try:
        lo, hi = (int(x) for x in spec.split("-"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected LO-HI, got {spec!r}")
    if lo > hi:
        raise argparse.ArgumentTypeError(f"LO > HI in {spec!r}")
    return lo, hi


def add_shared(sp, outputs):
    sp.add_argument("--bam", required=True, help="transcriptome alignments (reference names = transcript ids)")
    sp.add_argument("--gtf", required=True)
    sp.add_argument("--fasta", required=True, help="genome FASTA of the GTF")
    sp.add_argument("--read-lengths", required=True, type=read_lengths, metavar="LO-HI")
    sp.add_argument("--out-prefix", required=True, help=f"writes <prefix>.{outputs}")
    sp.add_argument("--offsets", help="use this offsets.tsv instead of estimating offsets from the BAM")
    sp.add_argument("--min-offset-support", type=float, default=30,
                    help="reads spanning a start or stop codon needed to give a read length offsets")
    sp.add_argument("--tol", type=float, default=1e-3,
                    help="EM stops when one step moves no expected count by this many reads")
    sp.add_argument("--max-iter", type=int, default=100_000,
                    help="give up after this many EM steps (SQUAREM counts 2-3 per cycle)")


def main(argv=None):
    p = argparse.ArgumentParser(prog="ribokit")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("quant", help="P-site offsets and CDS quantification from a transcriptome BAM")
    add_shared(q, "quant.tsv, .offsets.tsv, .ties.tsv, .stats.tsv, .psites.tsv")
    o = sub.add_parser("orfs", help="reads per ORF in its reading frame (one P-site offset per read length), "
                                    "and outside the ORFs per transcript")
    add_shared(o, "orfs.tsv, .frames.tsv, .codons.tsv, .offsets.tsv, .ties.tsv, .stats.tsv, "
                  ".psites.tsv (every P-site)")
    o.add_argument("--orfs", help="ORF table, TSV with columns ORF_id Name start end (0-based transcript "
                                  "coordinates: start = first nt of the start codon, end = one past the last "
                                  "sense codon). The annotated CDSs are added")
    s = sub.add_parser("score", help="frame evidence per ORF (codon lead), pooled over `ribokit orfs` runs")
    s.add_argument("--orfs-prefix", nargs="+", required=True, metavar="PREFIX",
                   help="out prefixes of `ribokit orfs` runs with the same ORF table; their codons are summed")
    s.add_argument("--out-prefix", required=True, help="writes <prefix>.scores.tsv, .decoys.tsv, .score_stats.tsv")
    s.add_argument("--frame-lengths", type=read_lengths, metavar="LO-HI",
                   help=f"read lengths that vote (default: those with a frame-0 share >= {score.MIN_FRAME0} "
                        "in every library)")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    if a.command == "quant":
        quantify(a.bam, a.gtf, a.fasta, a.read_lengths, a.out_prefix, offsets_path=a.offsets,
                 min_support=a.min_offset_support, tol=a.tol, max_iter=a.max_iter)
    elif a.command == "orfs":
        orfs.run(a.bam, a.gtf, a.fasta, a.read_lengths, a.out_prefix, orfs_path=a.orfs, offsets_path=a.offsets,
                 min_support=a.min_offset_support, tol=a.tol, max_iter=a.max_iter)
    else:
        score.run(a.orfs_prefix, a.out_prefix, frame_lengths=a.frame_lengths)


if __name__ == "__main__":
    main()
