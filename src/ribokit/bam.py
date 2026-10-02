"""Footprint alignments from a transcriptome BAM.

Kept: mapped, forward-strand (sense to the transcript), primary or secondary
alignments whose CIGAR is matches plus optional soft clips at the ends. A 5' soft
clip is taken as untemplated nucleotides added during library prep: the 5' end is
the first aligned base and the length is the read length minus the 5' clip.
Alignments with I/D/N/H/P are dropped and counted. Each read name is one read
(names are UMI-collapsed upstream), whatever its number of alignments.
"""
import logging
from dataclasses import dataclass

import numpy as np
import pysam

log = logging.getLogger(__name__)

_MATCH = {0, 7, 8}   # M, =, X
_SOFT = 4


@dataclass
class Alignments:
    read: np.ndarray     # read index (0..n_reads-1), one per alignment
    tx: np.ndarray       # annotation index of the transcript, -1 if it has no valid CDS
    pos5: np.ndarray     # 0-based transcript position of the 5' end
    length: np.ndarray
    n_reads: int

    def subset(self, mask):
        return Alignments(self.read[mask], self.tx[mask], self.pos5[mask], self.length[mask], self.n_reads)

    def count_reads(self):
        return int(np.unique(self.read).size)


def _footprint_length(cigar):
    """Read length minus the 5' soft clip, or None if the CIGAR is not matches + end soft clips."""
    first, last = 0, len(cigar)
    three = 0
    if cigar[0][0] == _SOFT:
        first = 1
    if last > first and cigar[-1][0] == _SOFT:
        three, last = cigar[-1][1], last - 1
    aligned = 0
    for op, n in cigar[first:last]:
        if op not in _MATCH:
            return None
        aligned += n
    return aligned + three if aligned else None


def read_bam(path, anno, stats):
    """All forward alignments with an acceptable CIGAR. stats gets the read counts."""
    index = anno.index()
    names = {}
    read, tx, pos5, length = [], [], [], []
    n_bad_cigar_aln = 0
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        ref_tx = np.array([index.get(r.split("|", 1)[0], -1) for r in bam.references], dtype=np.int64)
        ref_len = np.array(bam.lengths, dtype=np.int64)
        known = ref_tx >= 0
        bad = known & (ref_len != anno.tx_len[np.maximum(ref_tx, 0)])
        if bad.any():
            i = np.flatnonzero(bad)[0]
            raise ValueError(f"{bam.references[i]} is {ref_len[i]} nt in the BAM but "
                             f"{anno.tx_len[ref_tx[i]]} nt in the annotation; {bad.sum()} such transcripts")
        stats["cds_transcripts_in_bam_header"] = int(known.sum())
        ref_tx_list = ref_tx.tolist()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.is_reverse or a.is_supplementary:
                continue
            r = names.setdefault(a.query_name, len(names))
            fp = _footprint_length(a.cigartuples)
            if fp is None:
                n_bad_cigar_aln += 1
                continue
            read.append(r)
            tx.append(ref_tx_list[a.reference_id])
            pos5.append(a.reference_start)
            length.append(fp)
    stats["reads_forward"] = len(names)
    stats["alignments_dropped_cigar"] = n_bad_cigar_aln
    aln = Alignments(np.array(read, dtype=np.int64), np.array(tx, dtype=np.int64),
                     np.array(pos5, dtype=np.int64), np.array(length, dtype=np.int64), len(names))
    stats["reads_cigar_ok"] = aln.count_reads()
    log.info("%s: %d forward reads, %d with an acceptable alignment (%d alignments dropped for their CIGAR)",
             path, len(names), stats["reads_cigar_ok"], n_bad_cigar_aln)
    return aln
