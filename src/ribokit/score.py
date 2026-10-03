"""Frame evidence along an ORF: the codon lead score, from `ribokit orfs` outputs.

Votes. Every codon of an ORF but its start codon votes, from the P-sites of the score lengths in its
3 nt (codons.tsv, summed over libraries; each alignment counts 1). It leads when more of them are in
the ORF's frame than in either other frame. A peak is one codon, so one vote at most.

Null: the ORF is not translated. A codon's reads are then multinomial over its 3 positions, each
position weighted by its expected reads from the EM: the ORF's own density as background, which has
no frame, plus, inside the host CDS, the CDS's density x 3 pi_l(the position's frame in the CDS).
p_i = P(lead) under it, exact for the codon's read count. Across lengths and libraries the frame
probabilities are mixed by their reads in the codon. leads ~ Poisson-binomial(p_i):
z = (leads - sum p_i) / sqrt(sum p_i (1 - p_i)), p = its upper tail, q = BH over the ORFs scored.
min_p is the p if every codon with reads led: how far the ORF's length and reads let it go.

Decoys: each ORF shifted by +1 and +2 nt, scored the same way with the ORF's density as background.
A shifted copy is left out where it overlaps the host CDS in the CDS's frame.
"""
import logging

import numpy as np
import pandas as pd
from scipy import stats

from . import orfs as orfs_mod

log = logging.getLogger(__name__)

# score lengths: frame-0 share at least this in every library (B.2)
MIN_FRAME0 = 0.9
COLUMNS = ["codons", "codons_with_reads", "reads", "in_frame_share", "leads", "expected", "z", "p", "min_p"]
# per-ORF sums and their value for an ORF without votes (None: NA)
EMPTY = {"codons_with_reads": 0, "reads": 0, "in_frame": 0, "leads": 0, "expected": 0.0, "z": None, "p": 1.0,
         "min_p": 1.0}


def lead_probability(n, q):
    """P(c0 > max(c1, c2)) for n reads, multinomial over the frames with probabilities q."""
    k0 = np.arange(1, n + 1)
    m = n - k0
    r = q[1] / (q[1] + q[2]) if q[1] + q[2] > 0 else 0.0
    lo, hi = np.maximum(m - k0 + 1, 0), np.minimum(m, k0 - 1)   # c1 in [lo, hi] keeps c0 ahead of both
    ahead = np.where(hi >= lo, stats.binom.cdf(hi, m, r) - stats.binom.cdf(lo - 1, m, r), 0.0)
    return float((stats.binom.pmf(k0, n, q[0]) * ahead).sum())


def upper_tail(p, k):
    """P(X >= k), X the sum of independent Bernoulli(p)."""
    dist = np.ones(1)
    for x in p:
        dist = np.r_[dist * (1 - x), 0.0] + np.r_[0.0, dist * x]
    return min(float(dist[k:].sum()), 1.0)


def score_lengths(frames):
    """Lengths whose frame-0 share is >= MIN_FRAME0 in every library."""
    ok = set.intersection(*(set(f.loc[f["frame0"] >= MIN_FRAME0, "length"]) for f in frames))
    return sorted(int(x) for x in ok)


def read_library(prefix):
    dtype = {"ORF_id": str, "Name": str}
    return (pd.read_csv(f"{prefix}.codons.tsv", sep="\t", dtype=dtype),
            pd.read_csv(f"{prefix}.orfs.tsv", sep="\t", dtype=dtype),
            pd.read_csv(f"{prefix}.frames.tsv", sep="\t"))


def score(libraries, lengths):
    """libraries: (codons, orfs, frames) per library, from the same ORF table. Returns the scores of
    the ORFs (all but the annotated CDSs), those of their decoys, and stats."""
    first = libraries[0][1]
    is_orf = first["type"].isin(orfs_mod.ORF_TYPES)
    for _, o, _ in libraries[1:]:
        if not o.loc[o["type"].isin(orfs_mod.ORF_TYPES), ["ORF_id", "Name", "start", "end"]].reset_index(
                drop=True).equals(first.loc[is_orf, ["ORF_id", "Name", "start", "end"]].reset_index(drop=True)):
            raise ValueError("the libraries' orfs.tsv list different ORFs")
    targets = first[is_orf & (first["type"] != "CDS")][["ORF_id", "Name", "type", "start", "end"]].reset_index(drop=True)
    cds = first[first["type"] == "CDS"].set_index("Name")
    s, e = targets["start"].to_numpy(), targets["end"].to_numpy()
    cs = targets["Name"].map(cds["start"]).fillna(-1).to_numpy(np.int64)
    ce = targets["Name"].map(cds["end"]).fillna(-1).to_numpy(np.int64)
    n_codons = (e - s) // 3
    index = pd.Series(np.arange(len(targets)), index=targets["ORF_id"])

    rows = []   # per library, shift and P-site position: (shift, ORF, codon, frame, reads, q0, q1, q2)
    for codons, o, frames in libraries:
        density = (o["NumReads"].fillna(0) / o["Length"]).set_axis(o["ORF_id"])
        d_orf = density[targets["ORF_id"]].to_numpy()
        d_cds = targets["Name"].map(cds["ORF_id"]).map(density).fillna(0).to_numpy()
        w = orfs_mod.frame_weights(frames)
        lo = int(frames["length"].iloc[0])
        c = codons[codons["length"].isin(lengths)]
        n = c[["frame0", "frame1", "frame2"]].to_numpy().ravel()
        t = np.repeat(index[c["ORF_id"]].to_numpy(), 3)[n > 0]
        rel = (3 * np.repeat(c["codon"].to_numpy(), 3) + np.tile(np.arange(3), len(c)))[n > 0]
        length = np.repeat(c["length"].to_numpy(), 3)[n > 0]
        n = n[n > 0]
        for shift in range(3):
            j, f = np.divmod(rel - shift, 3)
            m = (rel >= shift) & (j >= 1) & (j < n_codons[t])
            tm, jm = t[m], j[m]
            pos = (s[tm] + shift + 3 * jm)[:, None] + np.arange(3)
            in_cds = (pos >= cs[tm, None]) & (pos < ce[tm, None])
            expect = d_orf[tm, None] + d_cds[tm, None] * w[length[m, None] - lo, (pos - cs[tm, None]) % 3] * in_cds
            total = expect.sum(axis=1, keepdims=True)
            q = np.where(total > 0, expect / np.where(total > 0, total, 1), 1 / 3)
            rows.append(np.column_stack([np.full(m.sum(), shift), tm, jm, f[m], n[m], q]))
    r = pd.DataFrame(np.concatenate(rows) if rows else np.zeros((0, 8)),
                     columns=["shift", "orf", "codon", "frame", "n", "q0", "q1", "q2"])
    for col in ("shift", "orf", "codon", "frame"):
        r[col] = r[col].astype(np.int64)
    for f in range(3):
        r[f"c{f}"] = r["n"] * (r["frame"] == f)
        r[f"q{f}"] *= r["n"]                 # mixed by reads when summed
    per_codon = r.groupby(["shift", "orf", "codon"])[["c0", "c1", "c2", "q0", "q1", "q2"]].sum()
    c = per_codon[["c0", "c1", "c2"]].to_numpy()
    q = per_codon[["q0", "q1", "q2"]].to_numpy()
    reads = c.sum(axis=1)
    per_codon["reads"] = reads
    per_codon["lead"] = c[:, 0] > c[:, 1:].max(axis=1)
    per_codon["p_lead"] = [lead_probability(int(round(x)), y / x) for x, y in zip(reads, q)]

    def summary(g):
        p = g["p_lead"].to_numpy()
        leads = int(g["lead"].sum())
        var = (p * (1 - p)).sum()
        return pd.Series({"codons_with_reads": len(g), "reads": g["reads"].sum(), "in_frame": g["c0"].sum(),
                          "leads": leads, "expected": p.sum(),
                          "z": (leads - p.sum()) / np.sqrt(var) if var > 0 else np.nan,
                          "p": upper_tail(p, leads), "min_p": p.prod()})

    per_orf = (per_codon.groupby(["shift", "orf"]).apply(summary) if len(per_codon) else
               pd.DataFrame(columns=list(EMPTY), index=pd.MultiIndex.from_tuples([], names=["shift", "orf"])))

    out = []
    for shift in range(3):
        keep = np.ones(len(targets), dtype=bool)
        if shift:   # leave out a shifted copy that overlaps the host CDS in its frame
            keep = ~((cs >= 0) & (s + shift < ce) & (e + shift > cs) & ((s + shift - cs) % 3 == 0))
        tab = targets[keep].assign(start=s[keep] + shift, end=e[keep] + shift, shift=shift, codons=n_codons[keep] - 1)
        got = per_orf.reindex(pd.MultiIndex.from_arrays([np.full(keep.sum(), shift), np.flatnonzero(keep)]))
        got = got.astype(float).fillna({k: v for k, v in EMPTY.items() if v is not None})
        for col, default in EMPTY.items():
            tab[col] = got[col].to_numpy() if default is None else got[col].to_numpy().astype(type(default))
        with np.errstate(invalid="ignore"):
            tab["in_frame_share"] = tab.pop("in_frame") / tab["reads"]
        out.append(tab)
    scores = out[0].drop(columns="shift")[["ORF_id", "Name", "type", "start", "end"] + COLUMNS]
    scores["q"] = stats.false_discovery_control(scores["p"].to_numpy()) if len(scores) else []
    decoys = pd.concat(out[1:], ignore_index=True)[["ORF_id", "Name", "type", "start", "end", "shift"] + COLUMNS]
    st = {"libraries": len(libraries), "frame_lengths": ",".join(map(str, lengths)), "orfs_scored": len(scores),
          "decoys_scored": len(decoys), "decoys_left_out_cds_frame": 2 * len(targets) - len(decoys)}
    return scores, decoys, st


def run(prefixes, out_prefix, frame_lengths=None):
    libraries = [read_library(p) for p in prefixes]
    if frame_lengths is None:
        lengths = score_lengths([f for _, _, f in libraries])
        if not lengths:
            raise ValueError(f"no read length has a frame-0 share >= {MIN_FRAME0} in every library; "
                             "choose them with --frame-lengths")
    else:
        lengths = list(range(frame_lengths[0], frame_lengths[1] + 1))
    for _, _, f in libraries:
        if not set(lengths) <= set(f["length"]):
            raise ValueError(f"score lengths {lengths} are not all in the read length window of `ribokit orfs`")
    log.info("score lengths: %s", lengths)
    scores, decoys, st = score(libraries, lengths)
    log.info("%s", st)
    scores.to_csv(f"{out_prefix}.scores.tsv", sep="\t", index=False, na_rep="NA")
    decoys.to_csv(f"{out_prefix}.decoys.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(list(st.items()), columns=["stat", "value"]).to_csv(
        f"{out_prefix}.score_stats.tsv", sep="\t", index=False)
    return scores, decoys, st
