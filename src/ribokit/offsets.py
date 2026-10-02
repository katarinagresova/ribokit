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

Adjacent windows differ in one phase's offset, so how well the data pin a phase
is measured per phase: the best window against the best window that gives this
phase another offset ("rival"). z = score lead / sqrt(reads counted by exactly
one of the two), a McNemar statistic. On eIF4E 4h, every class whose offset
changed between the even and the odd transcripts of a library had z < 4.
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


def counted(ds, phase, a, b):
    """Reads whose P-site, at the window's offset for their phase, is in the CDS."""
    d = np.asarray(ds)[phase]
    return (a <= d) & (d <= b)


def phase_z(windows, phase, a, b, w):
    """z of the best window over each phase's rival; NaN if no window changes that phase."""
    best = windows[0]
    best_in = counted(best[2], phase, a, b)
    z = []
    for ph in range(3):
        rival = next((win for win in windows[1:] if win[2][ph] != best[2][ph]), None)
        if rival is None:
            z.append(np.nan)
            continue
        disagree = w[best_in != counted(rival[2], phase, a, b)].sum()
        z.append((best[0] - rival[0]) / np.sqrt(disagree) if disagree > 0 else 0.0)
    return z


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
        s = aln.length == length
        z = phase_z(windows, phase[s], a[s], b[s], w[s]) if ok else [np.nan] * 3
        for ph in range(3):
            rows.append((length, ph, windows[0][2][ph] if ok else np.nan, support[k0 + ph], z[ph], total[k0 + ph]))
    df = pd.DataFrame(rows, columns=["length", "phase", "offset", "support", "z", "reads"])
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
