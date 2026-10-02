"""P-site offsets per (read length, phase) by CDS inclusion (Ahmed et al. 2019).

offset = nt from the read's 5' end to the first nt of its P-site codon.
phase = (5' end - CDS start) mod 3, so a P-site lands on a codon start iff
(offset + phase) % 3 == 0. The score of an offset is the number of reads whose
P-site it puts in [CDS start, start of last sense codon]; only reads that span
the start or the stop codon change it, so they are a class's support.

Within one read length the phases differ by one nt of 5' trimming, so their
offsets are three consecutive values, one per phase. Each length gets the window
of three consecutive offsets with the highest summed score (the design record's
"window over offset - phase"). RiboStan's window runs over (offset, phase) rows
instead and cannot express e.g. (12, 11, 13); choosing each phase on its own
was unstable on real data (eIF4E 4h: 17 of 39 classes changed between libraries).
"""
import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# The P-site codon must leave at least one codon (E-site) before it and one
# (A-site) after it in the read.
MIN_FLANK = 3


def read_weights(aln):
    """1 / number of alignments of the read, so a multimapper counts once in total."""
    return 1.0 / np.bincount(aln.read, minlength=aln.n_reads)[aln.read]


def estimate_offsets(aln, anno, lengths, min_support):
    lo, hi = lengths
    t = aln.tx
    phase = (aln.pos5 - anno.cds_start[t]) % 3
    k = (aln.length - lo) * 3 + phase
    n_class, n_delta = (hi - lo + 1) * 3, hi + 1
    w = read_weights(aln)
    a = anno.cds_start[t] - aln.pos5         # smallest offset with the P-site in the CDS
    b = anno.cds_end[t] - 3 - aln.pos5       # largest
    first, last = np.maximum(a, 0), np.minimum(b, hi)
    ok = first <= last
    # score[k, d] = weighted reads of class k whose P-site at offset d is in the CDS,
    # from a difference array over d
    diff = (np.bincount(k[ok] * (n_delta + 1) + first[ok], weights=w[ok], minlength=n_class * (n_delta + 1))
            - np.bincount(k[ok] * (n_delta + 1) + last[ok] + 1, weights=w[ok], minlength=n_class * (n_delta + 1)))
    score = np.cumsum(diff.reshape(n_class, n_delta + 1), axis=1)[:, :n_delta]
    spans = ((a >= 0) & (a < aln.length)) | ((b >= 0) & (b < aln.length))
    support = np.bincount(k, weights=w * spans, minlength=n_class)
    total = np.bincount(k, weights=w, minlength=n_class)

    rows = []
    for length in range(lo, hi + 1):
        k0 = (length - lo) * 3
        windows = []   # (summed score, centre, offset per phase), window inside the allowed offsets
        for c in range(MIN_FLANK + 1, length - 3 - MIN_FLANK):
            ds = [next(d for d in (c - 1, c, c + 1) if (d + ph) % 3 == 0) for ph in range(3)]
            windows.append((sum(score[k0 + ph, d] for ph, d in enumerate(ds)), c, ds))
        windows.sort(key=lambda win: (-win[0], win[1]))
        sup = support[k0:k0 + 3].sum()
        ok = sup >= min_support and len(windows) > 0
        margin = (windows[0][0] - windows[1][0]) / sup if ok and len(windows) > 1 else np.nan
        for ph in range(3):
            rows.append((length, ph, windows[0][2][ph] if ok else np.nan, support[k0 + ph], margin, total[k0 + ph]))
    df = pd.DataFrame(rows, columns=["length", "phase", "offset", "support", "margin", "reads"])
    df["offset"] = df["offset"].astype("Int64")
    log.info("offsets for %d of %d (length, phase) classes", df["offset"].notna().sum(), len(df))
    return df


def read_offsets(path):
    df = pd.read_csv(path, sep="\t")
    df["offset"] = df["offset"].astype("Int64")
    return df


def lookup(df, lengths):
    """Array indexed by (length - lo) * 3 + phase: the offset, or -1 if the class has none."""
    lo, hi = lengths
    table = np.full((hi - lo + 1) * 3, -1, dtype=np.int64)
    for length, ph, off in df[["length", "phase", "offset"]].itertuples(index=False):
        if lo <= length <= hi and not pd.isna(off):
            table[(length - lo) * 3 + ph] = off
    return table
