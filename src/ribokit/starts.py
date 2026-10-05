"""Start sites from harringtonine P-sites, against the P-sites of matched elongation libraries.

Harringtonine stops ribosomes at the start codon, and the ribosomes that started before it run off.
A used start has more harringtonine P-sites in its window (the start codon +-1 codon) than the
elongation profile of the same transcript expects. A pause has the same peak in both.

Layout. The transcripts sit end to end in one array, GAP empty nt apart, so that no window reaches
the next transcript. Each position has a region: leader, CDS, trailer (the stop codon and after) or
transcript (no annotated CDS), and a distance in the region: nt to the CDS start in the leader, codons
from the CDS start in the CDS, nt from the stop codon in the trailer, nt from the 5' end in a transcript.

Start lengths. Per read length, the P-sites around the annotated starts. The start lengths are the
longest run of lengths whose harringtonine P-sites peak at nt 0 of the start codon in each library.
Only their P-sites count, one per alignment.

Background, per pool. The harringtonine / elongation ratio depends on the region, and in the CDS on
the distance from the start (run-off). h = s_t f_b e: a factor per bin of region and distance (BINS),
fit on null windows, and a scale per transcript, over its background. A null window is the window of
a null codon: no start codon (START_CODONS) within 8 nt, outside the window of an annotated start and
its codons 1-10 (RAISED, never background), with mostly unique reads (processed pseudogenes mirror
their parent's start peak).

The test. A window has h harringtonine and e elongation reads. Its elongation rate has a gamma prior:
mean lam0 (the sum over its nt of the elongation density of their zone in the transcript: leader, CDS
frame 0, 1 or 2, trailer, transcript), shape a per region (from the null windows). Then h ~ NB(mean R max(e, (a + e) / (a / lam0 + 1)), shape a + e),
with R = f_b s_t: the prior lifts a window with few elongation reads, but does not shrink a pause.
Pools add their means and variances. An extra overdispersion phi per region and bin of
expected reads makes at most 1% of the null windows reach p <= 0.01. q is BH per class: the annotated
starts, and the other candidates per region.

Candidates: the starts of an ORF table and the annotated starts, or a scan: each start codon with
harringtonine reads in its window (BH counts every start codon of the class: the others have p = 1).
The codons 1-10 of an annotated start are not called. Then the scan calls starts in the order of p:
a start with q < 0.05 is called if no called start stops it, and it stops the candidates within 8 nt
and in its codons 1-10 (harringtonine raises them). Each called start gives an ORF to its first
in-frame stop, typed against the annotated CDS.
"""
import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import betainc, gammainc

from . import annotation

log = logging.getLogger(__name__)

GAP = 64
REGIONS = ("leader", "CDS", "trailer", "transcript")
CLASSES = REGIONS + ("annotated",)   # BH classes: the annotated starts, the other candidates by region
START_CODONS = ("ATG", "CTG", "GTG", "TTG", "AAG", "ACG", "AGG", "ATA", "ATC", "ATT")
WIN = (-1, 3)        # a start's window: nt from its first nt, the start codon and the nt before it
RAISED = (-3, 33)    # an annotated start's window and its codons 1-10: harringtonine raises them
NEAR = (-8, 8)       # a null codon has no start codon this near: their windows would overlap
PROFILE = (-15, 63)  # nt from the annotated start: the kernel to codon 10, then codons 11-20 scale the elongation
KERNEL_END = 33
STRIDE = WIN[1] - WIN[0]   # one null codon per window width for the null windows: they do not overlap
MAX_MULTI = 0.5      # null windows: multimapped share of the window's reads below this
PRIOR_READS = 10     # pseudo-reads: transcript scale at the global scale, bin factor at 1
PRIOR_NT = 100       # pseudo-nt at the region's mean elongation density, per transcript and region
DEPTH = (0, 0.3, 1, 3, 10, 30)               # lower edges of the bins of expected reads in the window
TAIL = 0.01                                  # phi: at most this share of null windows at p <= TAIL
PHI_GRID = np.r_[0, np.geomspace(1e-3, 100, 51)]
MIN_NULL = 200       # null windows for a phi of its own: else the region's, else all windows'
# lower edges of the distance bins per region (the region's distance unit)
BINS = {
    "leader": (0, 30, 100, 300),
    "CDS": (0, 20, 50, 100, 200, 400, 700, 1000),
    "trailer": (0, 100, 300, 1000),
    "transcript": (0,),
}


@dataclass
class Layout:
    names: np.ndarray       # transcript ids, sorted
    first: np.ndarray       # flat position of each transcript's nt 0
    tx_len: np.ndarray
    cds_start: np.ndarray   # -1 without an annotated CDS
    cds_end: np.ndarray     # one past the last sense codon; -1 without an annotated CDS
    seq: np.ndarray         # uint8 ASCII, N in the gaps
    size: int

    def positions(self):
        """Per flat position: transcript index, region index (REGIONS) and distance; -1 in the gaps."""
        tx = np.full(self.size, -1, dtype=np.int32)
        region = np.full(self.size, -1, dtype=np.int8)
        dist = np.full(self.size, -1, dtype=np.int32)
        for i, (f, n, cs, ce) in enumerate(zip(self.first, self.tx_len, self.cds_start, self.cds_end)):
            tx[f:f + n] = i
            x = np.arange(n)
            if cs < 0:
                region[f:f + n], dist[f:f + n] = 3, x
            else:
                region[f:f + n] = np.select([x < cs, x < ce], [0, 1], 2)
                dist[f:f + n] = np.select([x < cs, x < ce], [cs - x, (x - cs) // 3], x - ce)
        return tx, region, dist

    def cds_starts(self):
        """Flat positions of the annotated starts, and their transcript indices."""
        i = np.flatnonzero(self.cds_start >= 0)
        return self.first[i] + self.cds_start[i], i


ZONES = 6   # leader, CDS frame 0, 1 and 2, trailer, transcript


def zone_index(layout, tx, region):
    """Per flat position: leader 0, CDS frame 0-2 1-3, trailer 4, transcript 5; -1 in the gaps. The prior
    density of the elongation reads is per nt of a zone: a window holds one or two frame-0 nt of a CDS."""
    z = np.where(region < 0, -1, np.choose(np.maximum(region, 0), [0, 1, 4, 5])).astype(np.int8)
    cds = np.flatnonzero(region == 1)
    t = tx[cds]
    z[cds] = 1 + (cds - layout.first[t] - layout.cds_start[t]) % 3
    return z


def make_layout(gtf, fasta, names):
    """The layout of the transcripts `names` that are in the GTF."""
    anno = annotation.load_annotation(gtf, fasta)
    seqs = annotation.transcript_seqs(gtf, fasta, names)
    ids = np.array(sorted(seqs), dtype=object)
    idx = anno.index()
    c = np.array([idx.get(t, -1) for t in ids], dtype=np.int64)
    tx_len = np.array([len(seqs[t]) for t in ids], dtype=np.int64)
    first = GAP + np.r_[0, np.cumsum(tx_len + GAP)[:-1]]
    size = int(first[-1] + tx_len[-1] + GAP)
    seq = np.full(size, ord("N"), dtype=np.uint8)
    for t, f in zip(ids, first):
        seq[f:f + len(seqs[t])] = np.frombuffer(seqs[t].encode(), dtype=np.uint8)
    has = c >= 0
    log.info("layout: %d transcripts (%d with a CDS), %d nt", len(ids), has.sum(), size)
    return Layout(ids, first, tx_len, np.where(has, anno.cds_start[c], -1), np.where(has, anno.cds_end[c], -1),
                  seq, size)


def codon_mask(layout, codons=START_CODONS):
    """True at the positions where one of `codons` starts."""
    code = np.full(256, 4, dtype=np.int64)
    for i, b in enumerate(b"ACGT"):
        code[b] = i
    s = code[layout.seq]
    k = 25 * s[:-2] + 5 * s[1:-1] + s[2:]
    table = np.zeros(125, dtype=bool)
    for cod in codons:
        table[25 * "ACGT".index(cod[0]) + 5 * "ACGT".index(cod[1]) + "ACGT".index(cod[2])] = True
    return np.r_[table[k], False, False]


def covered(size, lo, hi):
    """True at the positions in any [lo, hi)."""
    d = np.zeros(size + 1, dtype=np.int64)
    np.add.at(d, lo, 1)
    np.add.at(d, hi, -1)
    return np.cumsum(d)[:-1] > 0


def background(layout, tx):
    """True at the transcript positions outside the RAISED zone of the annotated starts."""
    pos, _ = layout.cds_starts()
    return (tx >= 0) & ~covered(layout.size, pos + RAISED[0], pos + RAISED[1])


def null_codons(bg, starts, near=NEAR):
    """Positions c whose window is all background (bg) and that have no position of `starts` in
    [c + near[0], c + near[1]]. A gap is longer than a window, so the window is in c's transcript."""
    m, s = np.r_[0, np.cumsum(bg)], np.r_[0, np.cumsum(starts)]
    c = np.flatnonzero(bg)
    return c[(m[c + WIN[1]] - m[c + WIN[0]] == WIN[1] - WIN[0]) & (s[c + near[1] + 1] - s[c + near[0]] == 0)]


def read_psites(prefix, layout, lengths):
    """Flat positions, read lengths and multimapper flags (the read has another alignment) of the P-sites
    in <prefix>.psites.tsv (one per alignment) whose length is in [lo, hi] and whose transcript is in
    the layout."""
    first = pd.Series(layout.first, index=layout.names)
    ps = pd.read_csv(f"{prefix}.psites.tsv", sep="\t", dtype={"read": str, "Name": str, "psite": np.int64,
                                                              "length": np.int64})
    ps = ps[ps["length"].between(*lengths)]
    multi = ps["read"].duplicated(keep=False).to_numpy()
    f = ps["Name"].map(first).to_numpy()
    k = ~np.isnan(f)
    return (f[k].astype(np.int64) + ps["psite"].to_numpy()[k], ps["length"].to_numpy()[k].astype(np.int16),
            multi[k])


def window_sums(x, pos, win=WIN):
    """Sum of x over [pos + win[0], pos + win[1])."""
    c = np.r_[0, np.cumsum(x)]
    return c[pos + win[1]] - c[pos + win[0]]


def bin_index(region, dist):
    """Index of the (region, distance) bin, over BINS in order; -1 in the gaps."""
    out = np.full(region.shape, -1, dtype=np.int32)
    k = 0
    for r, name in enumerate(REGIONS):
        m = region == r
        out[m] = k + np.searchsorted(BINS[name], dist[m], "right") - 1
        k += len(BINS[name])
    return out


def bin_names():
    return [f"{name} {lo}" for name in REGIONS for lo in BINS[name]]


def fit_factors(h, e, tx, b, fit, scale, n_tx, n_bins, prior_reads, iters=500, tol=1e-9):
    """Poisson fit of h = s_t f_b e by alternating updates (h, e, tx, b per position): the factor per bin f
    over the positions `fit`, scaled to an e-weighted mean of 1 there, and the scale per transcript s over
    the positions `scale`. Priors: prior_reads pseudo-reads at the global scale for s; one pseudo-read at
    factor 1 for f, so a bin without reads does not get 0."""
    def sums(x, at):
        return np.bincount(tx[at].astype(np.int64) * n_bins + b[at], weights=x[at],
                           minlength=n_tx * n_bins).reshape(n_tx, n_bins)
    hs, es, hf, ef = sums(h, scale), sums(e, scale), sums(h, fit), sums(e, fit)
    f = np.ones(n_bins)
    for _ in range(iters):
        fe = (es * f).sum(axis=1)
        s = (hs.sum(axis=1) + prior_reads * hs.sum() / fe.sum()) / (fe + prior_reads)
        new = (hf.sum(axis=0) + 1) / ((s[:, None] * ef).sum(axis=0) + 1)
        new /= (new * ef.sum(axis=0)).sum() / ef.sum()
        done = np.abs(new - f).max() < tol
        f = new
        if done:
            break
    return f, s


def nb_sf(h, mean, size):
    """P(X >= h), X negative binomial with this mean and size (variance mean + mean^2 / size);
    Poisson if size is inf."""
    h = np.asarray(h, dtype=np.int64)
    mean, size = np.broadcast_to(mean, h.shape).astype(float), np.broadcast_to(size, h.shape).astype(float)
    out = np.ones(h.shape)
    k = h > 0
    poi = k & np.isinf(size)
    out[poi] = gammainc(h[poi], mean[poi])
    nb = k & ~np.isinf(size)
    out[nb] = betainc(h[nb], size[nb], mean[nb] / (size[nb] + mean[nb]))
    return out


def profile(pos, length, tx, cs_flat, lengths):
    """P-sites per (read length in [lo, hi], nt from the annotated start in PROFILE)."""
    lo, hi = lengths
    t = tx[pos]
    c = cs_flat[np.maximum(t, 0)]
    d = pos - c
    k = (t >= 0) & (c >= 0) & (d >= PROFILE[0]) & (d < PROFILE[1]) & (length >= lo) & (length <= hi)
    n = PROFILE[1] - PROFILE[0]
    return np.bincount((length[k] - lo) * n + d[k] - PROFILE[0], minlength=(hi - lo + 1) * n).reshape(-1, n)


def choose_start_lengths(peaks):
    """The longest run of consecutive lengths that peak at nt 0 in every library (peaks: length x library)."""
    best, cur = [], []
    for length, row in peaks.iterrows():
        cur = cur + [int(length)] if (row == 0).all() else []
        best = cur if len(cur) > len(best) else best
    return best


def kernel(hp, ep):
    """Per read length: harringtonine minus elongation P-sites (scaled on codons 11-20), as a share of
    nt -15..32, from profiles at the annotated starts."""
    far = slice(KERNEL_END - PROFILE[0], PROFILE[1] - PROFILE[0])
    k = hp[:, far].sum(axis=1) / np.maximum(ep[:, far].sum(axis=1), 1)
    kern = np.maximum(hp - k[:, None] * ep, 0)[:, :KERNEL_END - PROFILE[0]]
    return kern / np.maximum(kern.sum(axis=1, keepdims=True), 1e-300)


@dataclass
class Pool:
    h: np.ndarray           # harringtonine P-sites of the start lengths per position, summed over the pool
    e: np.ndarray           # elongation
    multi: np.ndarray       # multimapped P-sites, both kinds
    f: np.ndarray = None    # factor per bin
    s: np.ndarray = None    # scale per transcript
    dens: np.ndarray = None  # prior mean elongation reads per nt, per (transcript, zone)
    a: np.ndarray = None    # prior shape per region (inf: no spread)


def moments_nb(h, m, base_var):
    """Moment estimate of phi in var = base_var + phi m^2, at least 0."""
    return float(max((((h - m) ** 2).sum() - base_var.sum()) / max((m ** 2).sum(), 1e-300), 0))


def fit_pool(pool, tx, region, zone, bins, bg, null_pos, null_win, n_tx):
    """The pool's background: factors and scales (fit_factors), the prior mean and shape of the
    elongation reads of a window."""
    bg_pos = np.flatnonzero(bg)
    pool.f, pool.s = fit_factors(pool.h, pool.e, tx, bins, null_pos, bg_pos, n_tx, len(bin_names()), PRIOR_READS)
    key = tx[bg_pos].astype(np.int64) * ZONES + zone[bg_pos]
    z_mean = (np.bincount(zone[bg_pos], weights=pool.e[bg_pos], minlength=ZONES)
              / np.maximum(np.bincount(zone[bg_pos], minlength=ZONES), 1))
    pool.dens = ((np.bincount(key, weights=pool.e[bg_pos], minlength=n_tx * ZONES) + PRIOR_NT * np.tile(z_mean, n_tx))
                 / (np.bincount(key, minlength=n_tx * ZONES) + PRIOR_NT))
    e = window_sums(pool.e, null_win)
    lam0 = prior_mean(pool, null_win, tx, zone)
    pool.a = np.full(4, np.inf)
    for r in range(4):
        k = region[null_win] == r
        inv = moments_nb(e[k].astype(float), lam0[k], lam0[k]) if k.any() else 0.0
        pool.a[r] = 1 / inv if inv > 0 else np.inf


def prior_mean(pool, pos, tx, zone):
    """Prior mean elongation reads of the windows at pos: the densities of their nt's zones in their transcript."""
    t = tx[pos].astype(np.int64)
    lam = np.zeros(len(pos))
    for k in range(WIN[0], WIN[1]):
        z = zone[pos + k]
        lam += np.where((z >= 0) & (tx[pos + k] == t), pool.dens[t * ZONES + np.maximum(z, 0)], 0.0)
    return lam


def window_moments(pool, pos, tx, region, zone, bins):
    """Harringtonine reads in the windows at pos, and their mean and variance without phi: the
    elongation rate of a window has a gamma prior (mean lam0, shape a), updated by its e reads. The
    mean is R times the larger of e and the posterior mean: the prior lifts a window with few elongation
    reads, but does not shrink an elongation peak (a pause)."""
    t, r = tx[pos].astype(np.int64), region[pos]
    e = window_sums(pool.e, pos)
    R = pool.f[bins[pos]] * pool.s[t]
    lam0 = prior_mean(pool, pos, tx, zone)
    a = pool.a[r]
    fin = np.isfinite(a)
    af = np.where(fin, a, 0.0)
    m = R * np.maximum(e, np.where(fin, (af + e) / (af / lam0 + 1), lam0))
    var = m + np.where(fin, m ** 2 / np.where(fin, af + e, 1), 0.0)
    return window_sums(pool.h, pos), e, m, var


def nb_size(m, var, phi):
    """NB size with mean m and variance var + phi m^2; inf (Poisson) if that is not above m."""
    extra = var - m + phi * m ** 2
    return np.where(extra > 0, m ** 2 / np.maximum(extra, 1e-300), np.inf)


def tail_phi(h, m, var):
    """The smallest phi on PHI_GRID with at most TAIL of the windows at p <= TAIL."""
    for phi in PHI_GRID:
        if (nb_sf(h, m, nb_size(m, var, phi)) <= TAIL).mean() <= TAIL:
            return float(phi)
    return float(PHI_GRID[-1])


def depth_bin(m):
    return np.searchsorted(DEPTH, m, "right") - 1


def calibrate_phi(h, m, var, region):
    """phi per (region, depth bin) from null windows; MIN_NULL or more windows each, else the region's,
    else all windows'."""
    d = depth_bin(m)
    every = tail_phi(h, m, var)
    phi = {}
    for r in range(4):
        k = region == r
        own = tail_phi(h[k], m[k], var[k]) if k.sum() >= MIN_NULL else every
        for b in range(len(DEPTH)):
            kb = k & (d == b)
            phi[r, b] = tail_phi(h[kb], m[kb], var[kb]) if kb.sum() >= MIN_NULL else own
    return phi


def window_phi(m, region, phi):
    return np.array([phi[r, b] for r, b in zip(region.tolist(), depth_bin(m).tolist())], dtype=float)


def pvalues(h, m, var, region, phi):
    """p of each window: h ~ NB(m, var + phi m^2) with the phi of its region and depth bin."""
    return nb_sf(h, m, nb_size(m, var, window_phi(m, region, phi)))


def given_starts(path, layout):
    """The starts of an ORF table (TSV ORF_id Name start), plus the annotated starts: flat position ->
    ORF_ids. Rows whose transcript is not in the layout or whose start codon is off it are counted."""
    table = pd.read_csv(path, sep="\t", dtype={"ORF_id": str, "Name": str})[["ORF_id", "Name", "start"]]
    idx = pd.Series(np.arange(len(layout.names)), index=layout.names)
    t = table["Name"].map(idx)
    t, s = t.fillna(0).to_numpy(np.int64), table["start"].to_numpy(np.int64)
    ok = table["Name"].isin(idx.index).to_numpy() & (s >= 0) & (s + 3 <= layout.tx_len[t])
    c = np.flatnonzero(layout.cds_start >= 0)
    ids = pd.concat([pd.Series(table["ORF_id"].to_numpy()[ok], index=layout.first[t[ok]] + s[ok]),
                     pd.Series([f"{n}:CDS" for n in layout.names[c]], index=layout.first[c] + layout.cds_start[c])])
    ids = ids.groupby(level=0, sort=True).agg(",".join)
    return ids.index.to_numpy(np.int64), ids.to_numpy(), {"orf_table_rows": len(table),
                                                          "orf_table_rows_dropped": int((~ok).sum())}


def scan_candidates(layout, tx, pools, cs_pos):
    """Flat positions of the start codons (START_CODONS) with harringtonine reads in their window, and the
    annotated starts."""
    c = np.flatnonzero(codon_mask(layout) & (tx >= 0))
    h = sum(window_sums(pl.h, c) for pl in pools.values())
    return np.union1d(c[h > 0], cs_pos)


def next_stops(layout):
    """Per flat position x: the first in-frame stop codon at x or after (same frame as x), or the array size."""
    stop = codon_mask(layout, tuple(sorted(annotation.STOPS)))
    out = np.empty(layout.size, dtype=np.int64)
    for r in range(3):
        i = np.arange(r, layout.size, 3)
        out[i] = np.minimum.accumulate(np.where(stop[i], i, layout.size)[::-1])[::-1]
    return out


def orf_types(start, end, cs, ce):
    """ORF type from transcript coordinates against the annotated CDS [cs, ce) (cs -1: no CDS)."""
    up, inside = start < cs, (start > cs) & (start < ce)
    frame0 = (start - cs) % 3 == 0
    return np.select([cs < 0, start == cs, up & (end <= cs), up & frame0, up, inside & frame0, inside],
                     ["other", "CDS", "uORF", "extension", "uoORF", "truncation", "internal"], "dORF")


def bh(p, m):
    """Benjamini-Hochberg q of p among m tests; the m - len(p) tests not given have p = 1."""
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * m / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1.0)
    return out


def near_starts(pos, p, passing, annotated=()):
    """Called starts. The candidates in codons 1-10 of an annotated start (indices `annotated`) are not called:
    harringtonine raises them when it is used. Then by p (then position), a passing candidate that no called start
    stops is called. It stops the candidates within NEAR nt and in its codons 1-10. Returns called, and per
    candidate the index of the called or annotated start that stops it (-1: none)."""
    called = np.zeros(len(pos), dtype=bool)
    stopped_by = np.full(len(pos), -1, dtype=np.int64)
    for i in annotated:
        lo, hi = np.searchsorted(pos, [pos[i] + 1, pos[i] + RAISED[1]])
        stopped_by[lo:hi] = i
    for i in np.lexsort((pos, p)):
        if not passing[i] or stopped_by[i] >= 0:
            continue
        called[i] = True
        lo, hi = np.searchsorted(pos, [pos[i] + NEAR[0], pos[i] + RAISED[1]])
        near = np.arange(lo, hi)
        near = near[~called[near] & (stopped_by[near] < 0)]
        stopped_by[near] = i
    return called, stopped_by


def run(sheet_path, gtf, fasta, out_prefix, orfs_path=None, scan=False, start_lengths=None):
    """Start evidence from the `ribokit orfs` runs of a library sheet (TSV pool, harringtonine, elongation:
    out prefixes, one row per replicate), for the starts of an ORF table and the annotated starts, or for every
    start codon with harringtonine reads (scan: also the called starts and their ORFs)."""
    if (orfs_path is None) == (not scan):
        raise ValueError("give an ORF table or scan, not both")
    st = {}
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str)
    prefixes = list(sheet["harringtonine"]) + list(sheet["elongation"])
    names = None
    for p in prefixes:
        n = set(pd.read_csv(f"{p}.orfs.tsv", sep="\t", usecols=["Name"], dtype=str)["Name"])
        names = n if names is None else names & n
    layout = make_layout(gtf, fasta, names)
    tx, region, dist = layout.positions()
    bins = bin_index(region, dist)
    zone = zone_index(layout, tx, region)
    n_tx = len(layout.names)
    bg = background(layout, tx)
    nulls = null_codons(bg, codon_mask(layout))
    cs_pos, cs_tx = layout.cds_starts()
    cs_flat = np.full(n_tx, -1, dtype=np.int64)
    cs_flat[cs_tx] = cs_pos
    st.update(pools=sheet["pool"].nunique(), libraries=len(prefixes), transcripts=n_tx,
              transcripts_with_cds=len(cs_tx), null_codons=len(nulls))

    psites = {p: read_psites(p, layout, (0, np.iinfo(np.int16).max)) for p in prefixes}
    lo = min(int(x[1].min()) for x in psites.values())
    hi = max(int(x[1].max()) for x in psites.values())
    prof = {p: profile(pos, length, tx, cs_flat, (lo, hi)) for p, (pos, length, _) in psites.items()}
    peaks = pd.DataFrame({p: np.argmax(prof[p][:, :KERNEL_END - PROFILE[0]], axis=1) + PROFILE[0]
                          for p in sheet["harringtonine"]}, index=range(lo, hi + 1))
    if start_lengths is None:
        start_lengths = choose_start_lengths(peaks)
        if not start_lengths:
            raise ValueError("no read length peaks at nt 0 of the annotated starts in every harringtonine "
                             "library; choose them with --start-lengths")
    else:
        start_lengths = list(range(start_lengths[0], start_lengths[1] + 1))
    st["start_lengths"] = ",".join(map(str, start_lengths))
    log.info("start lengths: %s", start_lengths)

    kern_rows = []
    pools = {}
    for name, rows in sheet.groupby("pool", sort=True):
        hp = sum(prof[p] for p in rows["harringtonine"])
        ep = sum(prof[p] for p in rows["elongation"])
        nt = np.arange(PROFILE[0], KERNEL_END)
        for length, kr in zip(range(lo, hi + 1), kernel(hp, ep)):
            i = length - lo
            kern_rows.append(pd.DataFrame({"pool": name, "length": length, "nt": nt,
                                           "harringtonine": hp[i, :len(nt)], "elongation": ep[i, :len(nt)],
                                           "kernel": kr}))
        c = {kind: np.zeros(layout.size, dtype=np.int64) for kind in ("harringtonine", "elongation", "multi")}
        for kind in ("harringtonine", "elongation"):
            for p in rows[kind]:
                pos, length, multi = psites[p]
                k = np.isin(length, start_lengths)
                c[kind] += np.bincount(pos[k], minlength=layout.size)
                c["multi"] += np.bincount(pos[k & multi], minlength=layout.size)
        pools[name] = Pool(c["harringtonine"], c["elongation"], c["multi"])
        st[f"harringtonine_psites {name}"] = int(c["harringtonine"].sum())
        st[f"elongation_psites {name}"] = int(c["elongation"].sum())
    del psites

    reads = sum(window_sums(pl.h, nulls) + window_sums(pl.e, nulls) for pl in pools.values())
    multi = sum(window_sums(pl.multi, nulls) for pl in pools.values())
    null_pos = nulls[multi < MAX_MULTI * np.maximum(reads, 1)]
    null_win = null_pos[(null_pos - layout.first[tx[null_pos]]) % STRIDE == 0]
    st.update(null_codons_unique=len(null_pos), null_windows=len(null_win))
    factors = []
    for name, pl in pools.items():
        fit_pool(pl, tx, region, zone, bins, bg, null_pos, null_win, n_tx)
        factors.append(pd.DataFrame({"pool": name, "bin": bin_names(), "factor": pl.f}))
        st.update({f"prior_shape {name} {REGIONS[r]}": pl.a[r] for r in range(4)})

    def pooled(pos):
        got = [window_moments(pl, pos, tx, region, zone, bins) for pl in pools.values()]
        return tuple(sum(g[i] for g in got) for i in range(4))

    nh, _, nm, nv = pooled(null_win)
    phi = calibrate_phi(nh, nm, nv, region[null_win])
    null_p = pvalues(nh, nm, nv, region[null_win], phi)
    st.update({f"phi {REGIONS[r]} {DEPTH[b]}": v for (r, b), v in phi.items()})

    if scan:
        pos = scan_candidates(layout, tx, pools, cs_pos)
        st["candidates"] = len(pos)
    else:
        pos, ids, got = given_starts(orfs_path, layout)
        st.update(got)
    h, e, m, v = pooled(pos)
    p = pvalues(h, m, v, region[pos], phi)
    ann = np.isin(pos, cs_pos) & (m >= 3)
    enrich = np.median(h[ann] / m[ann]) if ann.any() else np.nan
    st["typical_enrichment"] = enrich
    typical_p = (pvalues(np.ceil(enrich * m).astype(np.int64), m, v, region[pos], phi) if ann.any()
                 else np.full(len(pos), np.nan))
    t = tx[pos].astype(np.int64)
    nt = pos - layout.first[t]
    out = pd.DataFrame({
        "Name": layout.names[t], "start": nt,
        "codon": [layout.seq[x:x + 3].tobytes().decode() for x in pos],
        "region": np.array(REGIONS)[region[pos]],
        "frame": pd.Series((nt - layout.cds_start[t]) % 3, dtype="Int64").mask(layout.cds_start[t] < 0),
        "ORF_id": [f"{n}_{x + 1}" for n, x in zip(layout.names[t], nt)] if scan else ids, "harringtonine": h, "elongation": e, "expected": m,
        "enrichment": h / np.maximum(m, 1e-300), "p": p, "typical_p": typical_p,
        "multimapped": sum(window_sums(pl.multi, pos) for pl in pools.values())
        / np.maximum(sum(window_sums(pl.h, pos) + window_sums(pl.e, pos) for pl in pools.values()), 1),
    })
    # BH per class, over every candidate of the class: a scan tests only the start codons with harringtonine
    # reads, but the others are tests too (p = 1)
    cls = np.where(np.isin(pos, cs_pos), 4, region[pos])
    if scan:
        every = np.flatnonzero(codon_mask(layout) & (tx >= 0))
        every = every[~np.isin(every, cs_pos)]
        n_tests = np.r_[np.bincount(region[every], minlength=4), len(cs_pos)]
    else:
        n_tests = np.bincount(cls, minlength=5)
    out["q"] = np.nan
    for c in range(5):
        k = cls == c
        if k.any():
            q = bh(p[k], n_tests[c])
            out.loc[k, "q"] = q
            cut = p[k][q < 0.05].max() if (q < 0.05).any() else 0.0
            nr = region[null_win] == min(c, 1)            # annotated starts: the CDS null windows
            rate = (null_p[nr] <= cut).mean() if nr.any() and cut > 0 else 0.0
            st.update({f"tests {CLASSES[c]}": int(n_tests[c]), f"candidates {CLASSES[c]}": int(k.sum()),
                       f"q<0.05 {CLASSES[c]}": int((q < 0.05).sum()), f"null_rate_at_cut {CLASSES[c]}": rate,
                       f"expected_false {CLASSES[c]}": rate * n_tests[c]})
    cols = ["Name", "start", "codon", "region", "frame", "ORF_id", "harringtonine", "elongation", "expected",
            "enrichment", "p", "typical_p", "q", "multimapped"]
    if scan:
        called, stopped_by = near_starts(pos, p, out["q"].to_numpy() < 0.05, np.flatnonzero(cls == 4))
        stop = next_stops(layout)[pos + 3]
        has_stop = (stop < layout.size) & (tx[np.minimum(stop, layout.size - 1)] == tx[pos])
        end = np.where(has_stop, stop - layout.first[t], -1)
        out["called"] = called
        out["stopped_by"] = pd.Series(out["ORF_id"].to_numpy()[np.maximum(stopped_by, 0)]).where(stopped_by >= 0)
        out["end"] = pd.Series(end, dtype="Int64").mask(~has_stop)
        out["type"] = np.where(has_stop, orf_types(nt, end, layout.cds_start[t], layout.cds_end[t]), "no stop")
        cols += ["called", "stopped_by", "end", "type"]
        orfs_out = out[called & has_stop & (out["type"] != "CDS").to_numpy()]
        orfs_out = orfs_out.rename(columns={"codon": "start_codon"})[
            ["ORF_id", "Name", "start", "end", "type", "start_codon", "harringtonine", "elongation", "expected",
             "enrichment", "p", "q"]]
        st.update(called=int(called.sum()), called_without_stop=int((called & ~has_stop).sum()),
                  stopped=int((stopped_by >= 0).sum()))
        st.update({f"start_orfs {k}": int(v) for k, v in orfs_out["type"].value_counts().sort_index().items()})
        orfs_out.to_csv(f"{out_prefix}.start_orfs.tsv", sep="\t", index=False, na_rep="NA")
    out = out[cols]
    log.info("%s", st)
    out.to_csv(f"{out_prefix}.starts.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(kern_rows, ignore_index=True).to_csv(f"{out_prefix}.kernel.tsv", sep="\t", index=False)
    pd.concat(factors, ignore_index=True).to_csv(f"{out_prefix}.factors.tsv", sep="\t", index=False)
    pd.DataFrame(list(st.items()), columns=["stat", "value"]).to_csv(f"{out_prefix}.starts_stats.tsv", sep="\t",
                                                                     index=False)
    return out, st
