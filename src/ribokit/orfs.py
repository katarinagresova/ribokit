"""Reads per ORF in its reading frame: uORFs, uoORFs and other ORFs besides the annotated CDS.

P-sites. quant's offsets depend on the phase of the 5' end relative to the CDS, so every
P-site lands on a codon of the CDS: right for the CDS, but it erases any other
reading frame. Here each read length has one offset, its phase-0 offset (phase 0
is the untrimmed phase, a multiple of 3), and P-site = 5' end + offset on every
transcript, wherever it lands. A read 1-2 nt short at the 5' end then puts its
P-site in frame 1 or 2 of the CDS. How often that happens is the length's frame
profile pi_l, measured on reads well inside CDSs (frames.tsv).

Components. The ORFs are the ORF table's plus the annotated CDSs. Every position of
a transcript that is in no ORF belongs to one outside component: the leader (5' of
the annotated CDS), the trailer (3' of it, its stop codon included), or the whole
transcript if it has no CDS. So every read with a P-site is counted somewhere.
Alignments to BAM references that are not transcripts of the GTF are dropped: a
read tied between such a reference and a GTF transcript counts for the latter.

EM. A read is compatible with each component that holds one of its P-sites.
P(read | ORF k) = 3 pi_l(phi) / L_k with phi = (P-site - start_k) mod 3, and
P(read | outside component) = 1 / L. When all of a read's components are ORFs that
put it in the same frame, the frame term cancels and this is quant's model.
"""
import logging

import numpy as np
import pandas as pd
import pysam

from . import annotation, bam, offsets, quant

log = logging.getLogger(__name__)

# frame profiles use reads at least this many nt inside the CDS at both ends:
# no start or stop peaks
INTERIOR = 15
ORF_TYPES = ("CDS", "uORF", "uoORF", "other")
OUTSIDE_TYPES = ("leader", "trailer", "transcript")


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


def frame_weights(frames):
    """3 pi_l(f), indexed (length - lo, frame): the EM's frame term. One pseudo-read per frame,
    so no weight is 0 and a length without interior reads gets 1 (no frame information)."""
    reads = frames["reads"].to_numpy()[:, None]
    n = np.nan_to_num(frames[["frame0", "frame1", "frame2"]].to_numpy()) * reads
    return 3 * (n + 1) / (reads + 3)


def orf_table(path, gtf, fasta, anno, bam_len, stats):
    """The ORFs: the table's (TSV ORF_id Name start end, in transcript coordinates: start = first
    nt of the start codon, end = one past the last sense codon) and the annotated CDSs it does not
    list. A row is dropped and counted if its transcript is not in the GTF, its length is not a
    positive multiple of 3, it runs off the transcript or no stop codon follows. The start codon
    is not checked. The type comes from the coordinates against the annotated CDS: uORF ends
    before it, uoORF starts before it out of frame and ends in or after it."""
    idx = anno.index()
    cds = pd.DataFrame({"ORF_id": [f"{t}:CDS" for t in anno.tx], "Name": anno.tx, "start": anno.cds_start,
                        "end": anno.cds_end, "start_codon": anno.start_codon})
    if path is None:
        table = cds.iloc[:0]
    else:
        table = pd.read_csv(path, sep="\t", dtype={"ORF_id": str, "Name": str})[["ORF_id", "Name", "start", "end"]]
        seqs = annotation.transcript_seqs(gtf, fasta, set(table["Name"]))
        for name, seq in seqs.items():
            if name in bam_len and bam_len[name] != len(seq):
                raise ValueError(f"{name} is {bam_len[name]} nt in the BAM but {len(seq)} nt in the annotation")
        dropped = dict.fromkeys(("transcript_not_in_gtf", "not_multiple_of_3", "off_transcript", "no_stop"), 0)
        keep, start_codon = [], []
        for name, s, e in zip(table["Name"], table["start"], table["end"]):
            seq = seqs.get(name)
            if seq is None:
                why = "transcript_not_in_gtf"
            elif e <= s or (e - s) % 3:
                why = "not_multiple_of_3"
            elif s < 0 or e + 3 > len(seq):
                why = "off_transcript"
            elif seq[e:e + 3] not in annotation.STOPS:
                why = "no_stop"
            else:
                keep.append(True)
                start_codon.append(seq[s:s + 3])
                continue
            dropped[why] += 1
            keep.append(False)
        stats["orf_table_rows"] = len(table)
        stats.update({f"orfs_dropped_{k}": v for k, v in dropped.items()})
        table = table[keep].assign(start_codon=start_codon)
    listed = set(zip(table["Name"], table["start"], table["end"]))
    added = cds[[k not in listed for k in zip(cds["Name"], cds["start"], cds["end"])]]
    stats["orfs_cds_added"] = len(added)
    orfs = pd.concat([table, added], ignore_index=True)
    if orfs["ORF_id"].duplicated().any():
        raise ValueError(f"ORF_id {orfs.loc[orfs['ORF_id'].duplicated(), 'ORF_id'].iloc[0]} is not unique")
    c = np.array([idx.get(n, -1) for n in orfs["Name"]], dtype=np.int64)
    has = c >= 0
    cs, ce = np.where(has, anno.cds_start[c], -1), np.where(has, anno.cds_end[c], -1)
    s, e = orfs["start"].to_numpy(np.int64), orfs["end"].to_numpy(np.int64)
    orfs["type"] = np.select([has & (s == cs) & (e == ce), has & (e <= cs), has & (s < cs) & ((cs - s) % 3 != 0)],
                             ORF_TYPES[:3], "other")
    stats.update({f"orfs_{t}": int((orfs["type"] == t).sum()) for t in ORF_TYPES})
    return orfs[["ORF_id", "Name", "type", "start", "end", "start_codon"]].astype({"start": np.int64, "end": np.int64})


def outside_components(orfs, anno, refs, ref_len):
    """Per BAM reference, its positions in no ORF: leader and trailer, or transcript if it has no
    annotated CDS. Components without a position are left out. Length = their positions."""
    idx = anno.index()
    spans = {}
    for name, s, e in zip(orfs["Name"], orfs["start"], orfs["end"]):
        spans.setdefault(name, []).append((s, e))
    rows = []
    for name, n in zip(refs, ref_len):
        c = idx.get(name)
        regions = ([("transcript", 0, n)] if c is None else
                   [("leader", 0, anno.cds_start[c]), ("trailer", anno.cds_end[c], n)])
        covered = np.zeros(n, dtype=bool)
        for s, e in spans.get(name, ()):
            covered[s:e] = True
        for kind, a, b in regions:
            free = int(b - a - covered[a:b].sum())
            if free > 0:
                rows.append((f"{name}:{kind}", name, kind, a, b, free))
    return pd.DataFrame(rows, columns=["ORF_id", "Name", "type", "start", "end", "Length"])


def orf_hits(ref, psite, orf_ref, start, end):
    """(alignment, ORF) pairs whose P-site is in the ORF: start <= P-site < end, same reference."""
    key = ref * 2**32 + psite
    order = np.argsort(key, kind="stable")
    lo = np.searchsorted(key[order], orf_ref * 2**32 + start)
    n = np.searchsorted(key[order], orf_ref * 2**32 + end) - lo
    i = order[np.repeat(lo - np.cumsum(n) + n, n) + np.arange(n.sum())]
    return i, np.repeat(np.arange(n.size), n)


def run(bam_path, gtf, fasta, lengths, out_prefix, orfs_path=None, offsets_path=None, min_support=30,
        tol=1e-3, max_iter=100_000):
    stats = {}
    anno = annotation.load_annotation(gtf, fasta)
    stats.update({f"annotation_{k}": v for k, v in anno.stats.items()})
    aln = bam.read_bam(bam_path, anno, stats)
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as f:
        ref_len = np.array(f.lengths, dtype=np.int64)
    # BAM references that are not GTF transcripts are not counted. Dropped before anything counts a
    # read's alignments (read_weights), so this is the run on a BAM without them.
    gtf_ids = set(annotation.read_gtf(gtf))
    in_gtf = np.array([r in gtf_ids for r in aln.refs], dtype=bool)
    keep = in_gtf[aln.ref]
    stats.update(refs_not_in_gtf=int((~in_gtf).sum()), alignments_dropped_ref_not_in_gtf=int((~keep).sum()))
    aln = aln.subset(keep)
    stats["reads_dropped_ref_not_in_gtf"] = stats["reads_cigar_ok"] - aln.count_reads()
    log.info("%d BAM references are not in the GTF: %d alignments dropped, %d reads with no other",
             stats["refs_not_in_gtf"], stats["alignments_dropped_ref_not_in_gtf"],
             stats["reads_dropped_ref_not_in_gtf"])
    orfs = orf_table(orfs_path, gtf, fasta, anno, dict(zip(aln.refs, ref_len.tolist())), stats)

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
    on = psite < ref_len[aln.ref]   # a long 3' soft clip at the transcript's end can push it past
    stats["alignments_psite_off_transcript"] = int((~on).sum())
    aln, psite = aln.subset(on), psite[on]
    stats["reads_with_psite"] = aln.count_reads()

    frames = frame_profile(aln, anno, psite, offsets.read_weights(aln), off, lengths)
    log.info("frame profiles:\n%s", frames.to_string(index=False))

    # components: the ORFs, then the outside components
    outside = outside_components(orfs, anno, [r for r, k in zip(aln.refs, in_gtf) if k], ref_len[in_gtf])
    comps = pd.concat([orfs.assign(Length=orfs["end"] - orfs["start"]), outside], ignore_index=True)
    n_orf = len(orfs)
    ref_index = {r: i for i, r in enumerate(aln.refs)}
    orf_ref = np.array([ref_index.get(n, -1) for n in orfs["Name"]], dtype=np.int64)
    start = orfs["start"].to_numpy()
    i, k = orf_hits(aln.ref, psite, orf_ref, start, orfs["end"].to_numpy())
    phi = (psite[i] - start[k]) % 3
    w_orf = frame_weights(frames)[aln.length[i] - lo, phi]
    # P-sites in no ORF: the outside component of their region on their reference
    out_index = {(n, t): n_orf + j for j, (n, t) in enumerate(zip(outside["Name"], outside["type"]))}
    ref_comp = np.array([[out_index.get((r, t), -1) for t in OUTSIDE_TYPES] for r in aln.refs],
                        dtype=np.int64).reshape(-1, 3)
    m = np.ones(aln.read.size, dtype=bool)
    m[i] = False
    region = np.where(aln.tx[m] < 0, 2, np.where(psite[m] < anno.cds_start[aln.tx[m]], 0, 1))
    comp_out = ref_comp[aln.ref[m], region]
    assert (comp_out >= 0).all()

    unique, classes, counts = quant.equivalence_classes(
        np.r_[aln.read[i], aln.read[m]], np.r_[k, comp_out], len(comps), np.r_[w_orf, np.ones(m.sum())])
    stats.update(reads_unique=int(unique.sum()), reads_multi=int(counts.sum()), equivalence_classes=len(counts))
    alpha, iterations, change, loglik = quant.em(unique, classes, counts, comps["Length"].to_numpy(float),
                                                 tol, max_iter)
    stats.update(em_iterations=iterations, em_last_max_change=change, em_tol=tol, em_loglik=loglik)
    log.info("EM: %d iterations, last max change %.3g reads, log-likelihood %.6f", iterations, change, loglik)
    if change >= tol:
        raise RuntimeError(f"EM did not converge in {max_iter} iterations (last max change {change:.3g} reads)")

    has_reads = (unique > 0) | (np.diff(classes.T.tocsr().indptr) > 0)
    comps["NumReads"] = np.where(has_reads, alpha, np.nan)
    ties = quant.tie_groups(unique, classes)
    stats.update(components=len(comps), components_with_reads=int(has_reads.sum()), tie_groups=len(ties),
                 components_in_ties=sum(len(g) for g in ties))
    stats.update({f"numreads_{t}": float(alpha[comps["type"].to_numpy() == t].sum()) for t in ORF_TYPES + OUTSIDE_TYPES})

    # for `ribokit score`: P-sites per codon, frame and length of each ORF but the annotated CDSs,
    # one count per alignment. Codon n is the stop codon: the +1 and +2 nt decoys reach into it.
    sel = np.flatnonzero(orfs["type"].to_numpy() != "CDS")
    ci, ck = orf_hits(aln.ref, psite, orf_ref[sel], start[sel], orfs["end"].to_numpy()[sel] + 3)
    rel = psite[ci] - start[sel][ck]
    codons = (pd.DataFrame({"orf": sel[ck], "codon": rel // 3, "length": aln.length[ci], "frame": rel % 3})
              .groupby(["orf", "codon", "length", "frame"]).size().unstack("frame")
              .reindex(columns=range(3)).fillna(0).astype(np.int64).reset_index())
    codons.insert(0, "ORF_id", orfs["ORF_id"].to_numpy()[codons.pop("orf")])
    codons.columns = ["ORF_id", "codon", "length", "frame0", "frame1", "frame2"]

    comps.sort_values(["Name", "start", "end", "ORF_id"], kind="stable").to_csv(
        f"{out_prefix}.orfs.tsv", sep="\t", index=False, na_rep="NA")
    offsets_df.to_csv(f"{out_prefix}.offsets.tsv", sep="\t", index=False, na_rep="NA")
    frames.to_csv(f"{out_prefix}.frames.tsv", sep="\t", index=False, na_rep="NA")
    codons.to_csv(f"{out_prefix}.codons.tsv", sep="\t", index=False)
    ids = comps["ORF_id"].to_numpy()
    pd.DataFrame([(ids[t], ids[g[0]]) for g in ties for t in g],
                 columns=["ORF_id", "tie_group"]).to_csv(f"{out_prefix}.ties.tsv", sep="\t", index=False)
    pd.DataFrame({
        "read": np.array(aln.read_names, dtype=object)[aln.read],
        "Name": np.array(aln.refs, dtype=object)[aln.ref],
        "psite": psite,
        "length": aln.length,
    }).to_csv(f"{out_prefix}.psites.tsv", sep="\t", index=False)
    pd.DataFrame(list(stats.items()), columns=["stat", "value"]).to_csv(
        f"{out_prefix}.stats.tsv", sep="\t", index=False)
    return comps, frames, offsets_df, stats
