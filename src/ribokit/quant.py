"""CDS quantification: assign footprints by P-site, then EM over equivalence classes.

A read is compatible with a transcript's CDS when one of its alignments to that
transcript puts the P-site inside the CDS. Model: a read comes from CDS t with
probability theta_t, and from each of its L_t positions equally, so
P(read) = sum over compatible t of theta_t / L_t. RiboStan's Stan model has no
1/L_t, which splits shared reads by read share instead of density.
EM starts from an even split of every read over its compatible CDSs (no
randomness) and runs until one EM step moves no expected count by `tol` reads.
SQUAREM extrapolates along pairs of EM steps: plain EM needs 10^5-10^6 steps
where near-identical CDSs (histone paralogs, UTR variants) split their reads
by a handful of distinguishing reads.
"""
import logging

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import annotation, bam, offsets

log = logging.getLogger(__name__)


def assign(aln, anno, table, lengths):
    """(has_offset, psite_in_cds) per alignment."""
    t = aln.tx
    phase = (aln.pos5 - anno.cds_start[t]) % 3
    off = table[(aln.length - lengths[0]) * 3 + phase]
    psite = aln.pos5 + off
    has_off = off >= 0
    return has_off, has_off & (psite >= anno.cds_start[t]) & (psite <= anno.cds_end[t] - 3)


def equivalence_classes(read, tx, n_orf):
    """unique[t] = reads compatible with t alone; classes = CSR (class x ORF), counts per class.
    A read with several compatible P-sites on one transcript counts once for it."""
    pairs = np.unique(read * n_orf + tx)
    r, t = pairs // n_orf, pairs % n_orf
    starts = np.flatnonzero(np.r_[True, r[1:] != r[:-1]])
    sizes = np.diff(np.r_[starts, r.size])
    single = sizes == 1
    unique = np.bincount(t[starts[single]], minlength=n_orf).astype(float)
    counts = {}
    for s, n in zip(starts[~single], sizes[~single]):
        key = t[s:s + n].tobytes()
        counts[key] = counts.get(key, 0) + 1
    members = [np.frombuffer(k, dtype=np.int64) for k in counts]
    indptr = np.r_[0, np.cumsum([m.size for m in members])].astype(np.int64)
    indices = np.concatenate(members) if members else np.zeros(0, dtype=np.int64)
    classes = sp.csr_matrix((np.ones(indices.size), indices, indptr), shape=(len(members), n_orf))
    return unique, classes, np.array(list(counts.values()), dtype=float)


def em(unique, classes, counts, lengths, tol, max_iter):
    """Expected reads per ORF, by EM accelerated with SQUAREM (Varadhan & Roland 2008, scheme S3).
    Stops when one EM step moves no expected count by `tol` reads; `max_iter` caps the EM steps.
    Returns (alpha, EM steps, last max change, log-likelihood)."""
    classes_t = classes.T.tocsr()
    sizes = np.diff(classes.indptr)
    n = unique.sum() + counts.sum()

    def step(alpha):
        w = alpha / lengths
        return unique + w * (classes_t @ (counts / (classes @ w)))

    def loglik(alpha):
        w = alpha / lengths / n
        return (unique[unique > 0] * np.log(w[unique > 0])).sum() + (counts * np.log(classes @ w)).sum()

    alpha = unique + classes_t @ (counts / sizes)
    ll, steps, change = loglik(alpha), 0, np.inf
    while steps < max_iter:
        a1 = step(alpha)
        steps += 1
        r = a1 - alpha
        change = np.max(np.abs(r)) if r.size else 0.0
        if change < tol:
            alpha = a1
            break
        a2 = step(a1)
        steps += 1
        v = a2 - 2 * a1 + alpha
        # sums, not BLAS dot products: their result can depend on the thread count
        vv = (v * v).sum()
        s = min(-np.sqrt((r * r).sum() / vv), -1.0) if vv > 0 else -1.0
        # s = -1 is a2, two plain EM steps; move s towards it until no positive count turns <= 0
        for _ in range(30):
            ext = alpha - 2 * s * r + s * s * v
            if (ext[alpha > 0] > 0).all():
                break
            s = (s - 1) / 2
        else:
            ext = a2
        a3 = step(ext)
        steps += 1
        ll3 = loglik(a3)
        if ll3 >= ll:
            alpha, ll = a3, ll3
        else:
            alpha, ll = a2, loglik(a2)
    return alpha, steps, change, loglik(alpha)


def tie_groups(unique, classes):
    """ORFs that no read tells apart: no unique reads and identical class membership."""
    classes_t = classes.T.tocsr()
    groups = {}
    for t in np.flatnonzero((unique == 0) & (np.diff(classes_t.indptr) > 0)):
        key = classes_t.indices[classes_t.indptr[t]:classes_t.indptr[t + 1]].tobytes()
        groups.setdefault(key, []).append(t)
    return [g for g in groups.values() if len(g) > 1]


def quantify(bam_path, gtf, fasta, lengths, out_prefix, offsets_path=None, min_support=30,
             tol=1e-3, max_iter=100_000):
    stats = {}
    anno = annotation.load_annotation(gtf, fasta)
    stats.update({f"annotation_{k}": v for k, v in anno.stats.items()})
    aln = bam.read_bam(bam_path, anno, stats)

    lo, hi = lengths
    aln = aln.subset((aln.length >= lo) & (aln.length <= hi))
    stats["reads_in_length_window"] = aln.count_reads()
    aln = aln.subset(aln.tx >= 0)
    stats["reads_on_cds_transcript"] = aln.count_reads()

    if offsets_path is None:
        offsets_df = offsets.estimate_offsets(aln, anno, lengths, min_support)
    else:
        offsets_df = offsets.read_offsets(offsets_path)
    table = offsets.lookup(offsets_df, lengths)
    has_off, in_cds = assign(aln, anno, table, lengths)
    stats["reads_with_offset"] = aln.subset(has_off).count_reads()
    aln = aln.subset(in_cds)
    stats["reads_assigned"] = aln.count_reads()

    unique, classes, counts = equivalence_classes(aln.read, aln.tx, len(anno.tx))
    stats.update(reads_unique=int(unique.sum()), reads_multi=int(counts.sum()), equivalence_classes=len(counts))
    alpha, iterations, change, loglik = em(unique, classes, counts, anno.cds_len.astype(float), tol, max_iter)
    stats.update(em_iterations=iterations, em_last_max_change=change, em_tol=tol, em_loglik=loglik)
    log.info("EM: %d iterations, last max change %.3g reads, log-likelihood %.6f", iterations, change, loglik)
    if change >= tol:
        raise RuntimeError(f"EM did not converge in {max_iter} iterations (last max change {change:.3g} reads)")

    has_reads = (unique > 0) | (np.diff(classes.T.tocsr().indptr) > 0)
    ties = tie_groups(unique, classes)
    stats.update(orfs=len(anno.tx), orfs_with_reads=int(has_reads.sum()), tie_groups=len(ties),
                 orfs_in_ties=sum(len(g) for g in ties))

    length = anno.cds_len
    density = alpha / length
    quant = pd.DataFrame({
        "Name": anno.tx,
        "Length": length,
        "EffectiveLength": length,
        "ritpm": np.where(has_reads, density / density.sum() * 1e6, np.nan),
        "NumReads": np.where(has_reads, alpha, np.nan),
    })
    quant.to_csv(f"{out_prefix}.quant.tsv", sep="\t", index=False, na_rep="NA")
    offsets_df.to_csv(f"{out_prefix}.offsets.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame([(anno.tx[t], anno.tx[g[0]]) for g in ties for t in g],
                 columns=["Name", "tie_group"]).to_csv(f"{out_prefix}.ties.tsv", sep="\t", index=False)
    pd.DataFrame(list(stats.items()), columns=["stat", "value"]).to_csv(
        f"{out_prefix}.stats.tsv", sep="\t", index=False)
    return quant, offsets_df, stats
