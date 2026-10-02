"""One CDS per transcript, in transcript coordinates, from a GTF and a genome FASTA.

Coordinates are 0-based. cds_start is the first nt of the start codon, cds_end is
one past the last nt of the last sense codon (the stop codon is never included,
whatever the GTF does). A CDS is kept if its length is a multiple of 3 and the
codon after it is a stop. Non-ATG starts are kept.
"""
import logging
import os
import re
import tempfile
from dataclasses import dataclass, field

import numpy as np
import pysam

log = logging.getLogger(__name__)

STOPS = {"TAA", "TAG", "TGA"}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ATTR = {key: re.compile(rf'{key} "([^"]+)"') for key in ("transcript_id", "gene_id")}


@dataclass
class Annotation:
    tx: list            # transcript ids with a valid CDS, sorted
    gene: list
    tx_len: np.ndarray
    cds_start: np.ndarray
    cds_end: np.ndarray
    start_codon: list
    stats: dict = field(default_factory=dict)

    @property
    def cds_len(self):
        return self.cds_end - self.cds_start

    def index(self):
        return {t: i for i, t in enumerate(self.tx)}


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def read_gtf(path):
    """{transcript_id: dict(chrom, strand, gene, exons, cds)} with 0-based half-open intervals."""
    txs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            tid = _ATTR["transcript_id"].search(f[8])
            if tid is None:
                continue
            gid = _ATTR["gene_id"].search(f[8])
            t = txs.setdefault(tid.group(1), dict(chrom=f[0], strand=f[6], gene=gid.group(1) if gid else "",
                                                  exons=[], cds=[]))
            if (t["chrom"], t["strand"]) != (f[0], f[6]):
                raise ValueError(f"transcript {tid.group(1)} has features on two chromosomes or strands")
            t["exons" if f[2] == "exon" else "cds"].append((int(f[3]) - 1, int(f[4])))
    return txs


def _tx_coord(exons, strand, g):
    """Transcript coordinate of genomic position g (exons in transcript order), or None."""
    done = 0
    for s, e in exons:
        if s <= g < e:
            return done + (g - s if strand == "+" else e - 1 - g)
        done += e - s
    return None


def _open_fasta(path):
    """pysam.FastaFile, indexing a copy in a temp dir if there is no .fai beside it."""
    if os.path.exists(path + ".fai"):
        return pysam.FastaFile(path)
    tmp = tempfile.mkdtemp(prefix="ribokit_fai_")
    fai = os.path.join(tmp, os.path.basename(path) + ".fai")
    log.info("no %s.fai, indexing into %s", path, tmp)
    pysam.faidx(path, "--fai-idx", fai)
    return pysam.FastaFile(path, filepath_index=fai)


def load_annotation(gtf, fasta):
    txs = read_gtf(gtf)
    fa = _open_fasta(fasta)
    chroms = set(fa.references)
    stats = dict(transcripts=len(txs))

    rows = []  # (tid, gene, tx_len, cds_start, cds_end_as_annotated, seq)
    n_no_cds = n_not_contiguous = 0
    for tid, t in txs.items():
        if not t["cds"]:
            n_no_cds += 1
            continue
        if t["chrom"] not in chroms:
            raise ValueError(f"{t['chrom']} (transcript {tid}) is not in {fasta}")
        exons = sorted(t["exons"], reverse=t["strand"] == "-")   # transcript order
        seq = "".join(fa.fetch(t["chrom"], s, e) for s, e in sorted(t["exons"])).upper()
        if t["strand"] == "-":
            seq = revcomp(seq)
        cs = min(s for s, _ in t["cds"])
        ce = max(e for _, e in t["cds"])
        first, last = (cs, ce - 1) if t["strand"] == "+" else (ce - 1, cs)
        start, end = _tx_coord(exons, t["strand"], first), _tx_coord(exons, t["strand"], last)
        cds_len = sum(e - s for s, e in t["cds"])
        if start is None or end is None or end - start + 1 != cds_len:
            n_not_contiguous += 1
            continue
        rows.append((tid, t["gene"], len(seq), start, end + 1, seq))
    stats.update(no_cds=n_no_cds, cds_not_contiguous_in_transcript=n_not_contiguous)

    # Does the GTF's CDS include the stop codon? Majority vote over CDSs whose length
    # is a multiple of 3, as RiboStan does.
    in3 = [r for r in rows if (r[4] - r[3]) % 3 == 0]
    last_is_stop = sum(r[5][r[4] - 3:r[4]] in STOPS for r in in3)
    next_is_stop = sum(r[5][r[4]:r[4] + 3] in STOPS for r in in3)
    if max(last_is_stop, next_is_stop) <= len(in3) / 2:
        raise ValueError(f"neither the last codon ({last_is_stop}) nor the next codon ({next_is_stop}) "
                         f"is a stop in most of {len(in3)} CDSs: wrong FASTA for this GTF?")
    includes_stop = last_is_stop > next_is_stop
    stats["gtf_cds_includes_stop"] = includes_stop

    kept = []
    n_not3 = n_nostop = 0
    for tid, gene, tx_len, start, end, seq in rows:
        if includes_stop:
            end -= 3
        if (end - start) % 3 != 0 or end <= start:
            n_not3 += 1
        elif seq[end:end + 3] not in STOPS:
            n_nostop += 1
        else:
            kept.append((tid, gene, tx_len, start, end, seq[start:start + 3]))
    kept.sort()
    n_non_atg = sum(k[5] != "ATG" for k in kept)
    stats.update(cds_not_multiple_of_3=n_not3, cds_without_stop=n_nostop, cds_kept=len(kept),
                 cds_kept_non_atg_start=n_non_atg)
    log.info("annotation: %s", stats)
    return Annotation(
        tx=[k[0] for k in kept],
        gene=[k[1] for k in kept],
        tx_len=np.array([k[2] for k in kept], dtype=np.int64),
        cds_start=np.array([k[3] for k in kept], dtype=np.int64),
        cds_end=np.array([k[4] for k in kept], dtype=np.int64),
        start_codon=[k[5] for k in kept],
        stats=stats,
    )
