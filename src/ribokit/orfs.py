"""P-sites that keep the reading frame, for counting ORFs that are not the annotated CDS.

quant's offsets depend on the phase of the 5' end relative to the CDS, so every
P-site lands on a codon of the CDS: right for the CDS, but it erases any other
reading frame. Here each read length has one offset, its phase-0 offset (phase 0
is the untrimmed phase, a multiple of 3), and P-site = 5' end + offset on every
transcript, wherever it lands. A read 1-2 nt short at the 5' end then puts its
P-site in frame 1 or 2 of the CDS. How often that happens is the length's frame
profile, measured on reads well inside CDSs (frames.tsv).
"""
import logging

import numpy as np
import pandas as pd

from . import annotation, bam, offsets

log = logging.getLogger(__name__)

# frame profiles use reads at least this many nt inside the CDS at both ends:
# no start or stop peaks
INTERIOR = 15


def psite_offsets(df, lengths):
    """Array indexed by length - lo: the length's phase-0 offset, or -1 if it has none."""
    return offsets.lookup(df, lengths)[0::3]


def frame_profile(aln, anno, psite, w, off, lengths):
    """Per read length: its offset and the share of P-sites in each frame of the CDS
    (frame = (P-site - CDS start) mod 3), weighted w, from reads at least INTERIOR nt
    inside a CDS at both ends."""
    lo, hi = lengths
    on = aln.tx >= 0
    t, pos5, length = aln.tx[on], aln.pos5[on], aln.length[on]
    cs, ce = anno.cds_start[t], anno.cds_end[t]
    inside = (pos5 >= cs + INTERIOR) & (pos5 + length <= ce - INTERIOR)
    k = (length - lo) * 3 + (psite[on] - cs) % 3
    n = np.bincount(k[inside], weights=w[on][inside], minlength=(hi - lo + 1) * 3).reshape(-1, 3)
    reads = n.sum(axis=1)
    with np.errstate(invalid="ignore"):
        share = n / reads[:, None]
    return pd.DataFrame({
        "length": np.arange(lo, hi + 1),
        "offset": pd.array(np.where(off >= 0, off, None), dtype="Int64"),
        "reads": reads,
        "frame0": share[:, 0],
        "frame1": share[:, 1],
        "frame2": share[:, 2],
    })


def run(bam_path, gtf, fasta, lengths, out_prefix, offsets_path=None, min_support=30):
    stats = {}
    anno = annotation.load_annotation(gtf, fasta)
    stats.update({f"annotation_{k}": v for k, v in anno.stats.items()})
    aln = bam.read_bam(bam_path, anno, stats)

    lo, hi = lengths
    aln = aln.subset((aln.length >= lo) & (aln.length <= hi))
    stats["reads_in_length_window"] = aln.count_reads()
    on_cds = aln.subset(aln.tx >= 0)
    stats["reads_on_cds_transcript"] = on_cds.count_reads()

    # estimated on the same alignments as in quant, so the two give the same table
    if offsets_path is None:
        offsets_df = offsets.estimate_offsets(on_cds, anno, lengths, min_support)
    else:
        offsets_df = offsets.read_offsets(offsets_path)
    off = psite_offsets(offsets_df, lengths)
    aln = aln.subset(off[aln.length - lo] >= 0)
    psite = aln.pos5 + off[aln.length - lo]
    stats["reads_with_psite"] = aln.count_reads()

    frames = frame_profile(aln, anno, psite, offsets.read_weights(aln), off, lengths)
    log.info("frame profiles:\n%s", frames.to_string(index=False))

    offsets_df.to_csv(f"{out_prefix}.offsets.tsv", sep="\t", index=False, na_rep="NA")
    frames.to_csv(f"{out_prefix}.frames.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame({
        "read": np.array(aln.read_names, dtype=object)[aln.read],
        "Name": np.array(aln.refs, dtype=object)[aln.ref],
        "psite": psite,
        "length": aln.length,
    }).to_csv(f"{out_prefix}.psites.tsv", sep="\t", index=False)
    pd.DataFrame(list(stats.items()), columns=["stat", "value"]).to_csv(
        f"{out_prefix}.stats.tsv", sep="\t", index=False)
    return frames, offsets_df, stats
