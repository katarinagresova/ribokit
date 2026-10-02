"""A small genome, GTF and transcriptome BAM with known truth, for the tests.

Genes are built in transcript orientation from exons of four kinds (5'UTR only;
UTR + ATG + codons; codons; codons + stop + 3'UTR), so every isoform that takes
exons in index order has a valid CDS. Exons are placed on one chromosome, in
reverse order for minus-strand genes.

Footprints: P-sites uniform over each CDS's codons plus start and stop peaks
(proportional to density), some UTR footprints, read length and offset from
PLANTED. Each read is aligned exactly to every transcript that contains its
sequence.
"""
import numpy as np
import pysam

STOPS = ("TAA", "TAG", "TGA")
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]
_COMP = str.maketrans("ACGT", "TGCA")

# read length -> [(offset, probability)]; phase = -offset mod 3, so one offset per (length, phase).
# Length 29 is the (12, 12, 15) case of the design record in RiboStan's convention.
PLANTED = {28: [(12, .7), (11, .1), (13, .2)],
           29: [(12, .6), (11, .15), (13, .25)],
           30: [(12, .5), (13, .3), (14, .2)]}
TRUE_OFFSETS = {(length, -d % 3): d for length, ds in PLANTED.items() for d, _ in ds}

# name, strand, exons (kind, utr nt, codons, utr nt), isoforms {name: exon indices}, CDS reads per isoform
GENES = [
    ("gA", "+", [("start", 60, 40, 0), ("mid", 0, 50, 0), ("end", 0, 60, 80)], {"A1": [0, 1, 2]}, {"A1": 3000}),
    ("gB", "+", [("start", 50, 30, 0), ("mid", 0, 40, 0), ("end", 0, 80, 60)],
     {"B1": [0, 1, 2], "B2": [0, 2]}, {"B1": 2000, "B2": 3000}),
    ("gC", "-", [("utr5", 70, 0, 0), ("utr5", 90, 0, 0), ("start", 30, 30, 0), ("end", 0, 100, 70)],
     {"C1": [0, 2, 3], "C2": [1, 2, 3]}, {"C1": 1500, "C2": 1500}),
    ("gD", "-", [("start", 40, 50, 0), ("mid", 0, 30, 0), ("mid", 0, 40, 0), ("end", 0, 50, 50)],
     {"D1": [0, 1, 2, 3]}, {"D1": 2500}),
    ("gN", "+", [("utr5", 500, 0, 0)], {"N1": [0]}, {}),
]


def revcomp(s):
    return s.translate(_COMP)[::-1]


def _exon_seq(rng, kind, utr5, codons, utr3):
    rand = lambda n: "".join(rng.choice(list("ACGT"), n))
    cds = "".join(rng.choice(SENSE, codons))
    if kind == "utr5":
        return rand(utr5)
    if kind == "start":
        return rand(utr5) + "ATG" + cds
    if kind == "mid":
        return cds
    return cds + rng.choice(STOPS) + rand(utr3)


def build(rng):
    """Genome sequence and transcripts {name: dict(gene, strand, seq, cds_start, cds_end, gexons)}."""
    genome, txs = [], {}
    pos = 0

    def place(seq):
        nonlocal pos
        genome.append("".join(rng.choice(list("ACGT"), 150)))
        pos += 150
        genome.append(seq)
        start, pos = pos, pos + len(seq)
        return start, pos

    for gene, strand, exon_specs, isoforms, _ in GENES:
        seqs = [_exon_seq(rng, *e) for e in exon_specs]
        order = range(len(seqs)) if strand == "+" else reversed(range(len(seqs)))
        gexons = {}
        for i in order:
            gexons[i] = place(seqs[i] if strand == "+" else revcomp(seqs[i]))
        for name, idx in isoforms.items():
            seq = "".join(seqs[i] for i in idx)
            tx = dict(gene=gene, strand=strand, seq=seq, gexons=[gexons[i] for i in idx], cds_start=None, cds_end=None)
            starts = [k for k, i in enumerate(idx) if exon_specs[i][0] == "start"]
            if starts:
                k = starts[0]
                cs = sum(len(seqs[i]) for i in idx[:k]) + exon_specs[idx[k]][1]
                n_codons = 1 + sum(exon_specs[i][2] for i in idx[k:])
                tx.update(cds_start=cs, cds_end=cs + 3 * n_codons)
                assert seq[cs:cs + 3] == "ATG" and seq[tx["cds_end"]:tx["cds_end"] + 3] in STOPS
            txs[name] = tx
    genome.append("".join(rng.choice(list("ACGT"), 150)))
    return "".join(genome), txs


def _genomic_segments(tx, start, end):
    """Transcript interval [start, end) as genomic (start, end) pieces, 0-based."""
    out, done = [], 0
    for gs, ge in tx["gexons"]:
        n = ge - gs
        lo, hi = max(start, done) - done, min(end, done + n) - done
        if lo < hi:
            out.append((gs + lo, gs + hi) if tx["strand"] == "+" else (ge - hi, ge - lo))
        done += n
    return out


def write_reference(directory, genome, txs, cds_includes_stop=False):
    fa = directory / "genome.fa"
    fa.write_text(">chr1\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    pysam.faidx(str(fa))
    gtf = directory / ("annotation_with_stop.gtf" if cds_includes_stop else "annotation.gtf")
    lines = []
    for name, tx in txs.items():
        attr = f'gene_id "{tx["gene"]}"; transcript_id "{name}";'
        gs, ge = min(s for s, _ in tx["gexons"]), max(e for _, e in tx["gexons"])
        lines.append(f"chr1\tsynth\ttranscript\t{gs + 1}\t{ge}\t.\t{tx['strand']}\t.\t{attr}")
        for s, e in tx["gexons"]:
            lines.append(f"chr1\tsynth\texon\t{s + 1}\t{e}\t.\t{tx['strand']}\t.\t{attr}")
        if tx["cds_start"] is not None:
            end = tx["cds_end"] + (3 if cds_includes_stop else 0)
            for s, e in _genomic_segments(tx, tx["cds_start"], end):
                lines.append(f"chr1\tsynth\tCDS\t{s + 1}\t{e}\t.\t{tx['strand']}\t.\t{attr}")
    gtf.write_text("\n".join(lines) + "\n")
    return fa, gtf


def simulate_reads(rng, txs, peak_codons=5, utr_density=0.3, long_reads=200, softclip_fraction=0.05,
                   junk_cigar_reads=7, unknown_ref_reads=11):
    """[(name, seq, cigar, [(tx, pos5)])] plus truth: CDS reads per transcript."""
    lengths = list(PLANTED)
    reads, truth = [], {}

    def footprint(tx, psite, length=None):
        length = length or rng.choice(lengths)
        ds, ps = zip(*PLANTED[length]) if length in PLANTED else ((12,), (1.0,))
        pos5 = psite - rng.choice(ds, p=ps)
        if pos5 < 0 or pos5 + length > len(tx["seq"]):
            return None
        return tx["seq"][pos5:pos5 + length]

    for gene, _, _, _, abundance in GENES:
        for name, n in abundance.items():
            tx = txs[name]
            cs, ce = tx["cds_start"], tx["cds_end"]
            psites = list(cs + 3 * rng.integers(0, (ce - cs) // 3, n))
            # start and stop peaks hold peak_codons codons' worth of dwell, so they scale
            # with density (reads per codon), as ribosomes pausing at the ends do
            peak = int(peak_codons * n / ((ce - cs) // 3))
            psites += [cs] * peak + [ce - 3] * peak
            seqs = [s for s in (footprint(tx, p) for p in psites) if s]
            truth[name] = len(seqs)
            utr = [p for p in range(len(tx["seq"])) if p < cs or p > ce - 3]
            seqs += [s for s in (footprint(tx, p) for p in rng.choice(utr, int(utr_density * len(utr)))) if s]
            seqs += [s for s in (footprint(tx, cs + 3 * rng.integers(0, (ce - cs) // 3), 35)
                                 for _ in range(long_reads // len(abundance))) if s]
            reads += seqs

    out = []
    for i, seq in enumerate(reads):
        hits = [(t, p) for t, tx in txs.items() for p in _find_all(tx["seq"], seq)]
        name = f"sim_{i}_x{rng.integers(1, 4)}"
        if rng.random() < softclip_fraction:
            out.append((name, rng.choice(list("ACGT")) + seq, f"1S{len(seq)}M", hits))
        else:
            out.append((name, seq, f"{len(seq)}M", hits))
    for j in range(junk_cigar_reads):
        name, (t, p) = f"junk_{j}_x1", (next(iter(txs)), 100)
        out.append((name, txs[t]["seq"][p:p + 27], "10M1D17M", [(t, p)]))
    for j in range(unknown_ref_reads):
        out.append((f"unknown_{j}_x1", "A" * 28, "28M", [("unannotated", 5)]))
    return out, truth


def _find_all(seq, sub):
    i = seq.find(sub)
    while i != -1:
        yield i
        i = seq.find(sub, i + 1)


def write_bam(path, txs, reads):
    refs = list(txs) + ["unannotated"]
    lens = [len(tx["seq"]) for tx in txs.values()] + [1000]
    header = {"HD": {"VN": "1.4", "SO": "unsorted"}, "SQ": [{"SN": r, "LN": n} for r, n in zip(refs, lens)]}
    ref_id = {r: i for i, r in enumerate(refs)}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, seq, cigar, hits in reads:
            for k, (t, p) in enumerate(hits):
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.query_sequence, a.cigarstring = name, seq, cigar
                a.flag = 0 if k == 0 else 256
                a.reference_id, a.reference_start, a.mapping_quality = ref_id[t], p, 255
                a.set_tag("NH", len(hits))
                out.write(a)


def make_dataset(directory, seed=1):
    rng = np.random.default_rng(seed)
    genome, txs = build(rng)
    fa, gtf = write_reference(directory, genome, txs)
    _, gtf_with_stop = write_reference(directory, genome, txs, cds_includes_stop=True)
    reads, truth = simulate_reads(rng, txs)
    bam_path = directory / "reads.bam"
    write_bam(bam_path, txs, reads)
    return dict(fasta=fa, gtf=gtf, gtf_with_stop=gtf_with_stop, bam=bam_path, txs=txs, truth=truth, reads=reads)
