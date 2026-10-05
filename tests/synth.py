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


# ORF dataset (make_orf_dataset): single-exon transcripts on the + strand, ORFs besides the CDSs.
# Reads per source: a CDS or ORF (P-sites on its codons, frames from PLANTED), or background
# (BACKGROUND reads per nt, uniform over the positions in no ORF).
ORF_READS = {"U1:CDS": 6000, "U1_uORF": 400, "U2_main": 2000, "U2_uoORF": 300, "U3:CDS": 2000, "U3_uoORF": 300,
             "U4:CDS": 2000, "U5:CDS": 2000, "U5_uoORF": 600, "N1_ORF": 300}
# reads with their P-site on the start codon: U1_CUG has no others, U4's CDS has a start peak whose
# reads 1-2 nt long at the 5' end put their P-site in U4_uoORF's leader part (untranslated)
PEAKS = {"U1_CUG": 150, "U4:CDS": 150}
BACKGROUND = 0.4


def build_orf_transcripts(rng):
    """Transcripts {name: dict(seq, cds_start, cds_end)} and ORFs {id: (name, start, end)}, CDSs included."""
    rand = lambda n: "".join(rng.choice(list("ACGT"), n))
    codons = lambda n: "".join(rng.choice(SENSE, n))
    stop = lambda: str(rng.choice(STOPS))
    txs, orfs = {}, {}
    # U1: a translated uORF, an untranslated one and a CUG ORF, then the CDS
    seq = ""
    for orf_id, gap, start, n in (("U1_uORF", 30, "ATG", 20), ("U1_cand", 24, "ATG", 15), ("U1_CUG", 30, "CTG", 10)):
        seq += rand(gap)
        orfs[orf_id] = ("U1", len(seq), len(seq) + 3 + 3 * n)
        seq += start + codons(n) + stop()
    seq += rand(40)
    orfs["U1:CDS"] = ("U1", len(seq), len(seq) + 450)
    seq += "ATG" + codons(149) + stop() + rand(120)
    txs["U1"] = seq
    # U2, U3, U4, U5: a uoORF in frame 1 / 2 / 1 / 1 of the CDS, from 60 nt (U5: 150 nt) upstream of
    # it to 90 nt into it. Its stop spans CDS codons 30 and 31.
    for name, f, up in (("U2", 1, 60), ("U3", 2, 60), ("U4", 1, 60), ("U5", 1, 150)):
        cs = up + 90
        cds = ["ATG"] + list(rng.choice(SENSE, 149))
        cds[30], cds[31] = ("CTA", "AAA") if f == 1 else ("CCT", "AAA")
        s, e = cs - up + f, cs + 90 + f
        txs[name] = rand(s) + "ATG" + rand(cs - s - 3) + "".join(cds) + stop() + rand(100)
        orfs[f"{name}_uoORF"] = (name, s, e)
        orfs["U2_main" if name == "U2" else f"{name}:CDS"] = (name, cs, cs + 450)
    # P1: no CDS, holds 240 nt of U1's CDS (pseudogene-like)
    cs1 = orfs["U1:CDS"][1]
    txs["P1"] = rand(80) + txs["U1"][cs1 + 150:cs1 + 390] + rand(80)
    # N1: no CDS, a translated ORF
    orfs["N1_ORF"] = ("N1", 100, 178)
    txs["N1"] = rand(100) + "ATG" + codons(25) + stop() + rand(150)
    for name, s, e in orfs.values():
        assert (e - s) % 3 == 0 and txs[name][e:e + 3] in STOPS
    cds = {"U1": "U1:CDS", "U2": "U2_main", "U3": "U3:CDS", "U4": "U4:CDS", "U5": "U5:CDS"}
    return {name: dict(seq=seq, cds_start=orfs[cds[name]][1] if name in cds else None,
                       cds_end=orfs[cds[name]][2] if name in cds else None) for name, seq in txs.items()}, orfs


def outside_id(txs, orfs, name, p):
    """The outside component holding position p of transcript `name`, None if p is in an ORF."""
    if any(n == name and s <= p < e for n, s, e in orfs.values()):
        return None
    tx = txs[name]
    if tx["cds_start"] is None:
        return f"{name}:transcript"
    return f"{name}:leader" if p < tx["cds_start"] else f"{name}:trailer"


def make_orf_dataset(directory, seed=2, background_under_orfs=False):
    """Reference, BAM and ORF table (with rows to drop, and U2's CDS listed as U2_main).
    truth[component] = reads drawn from it whose P-site at the phase-0 offset (12) is in it
    (ORF: its span; outside component: its positions); the rest are counted in truth["leaked"].
    background_under_orfs: background over every position, as for the score; truth then counts
    the background drawn under ORFs as leaked."""
    rng = np.random.default_rng(seed)
    txs, orfs = build_orf_transcripts(rng)
    genome, pos = [], 0
    for tx in txs.values():
        genome += ["".join(rng.choice(list("ACGT"), 150)), tx["seq"]]
        tx.update(gene=f"g_{len(genome)}", strand="+", gexons=[(pos + 150, pos + 150 + len(tx["seq"]))])
        pos += 150 + len(tx["seq"])
    fa, gtf = write_reference(directory, "".join(genome) + "A" * 150, txs)

    sources = []   # (component, transcript, P-site)
    for orf_id, (name, s, e) in orfs.items():
        sources += [(orf_id, name, p) for p in s + 3 * rng.integers(0, (e - s) // 3, ORF_READS.get(orf_id, 0))]
    for orf_id, n in PEAKS.items():
        sources += [(orf_id, orfs[orf_id][0], orfs[orf_id][1])] * n
    for name, tx in txs.items():
        free = [p for p in range(len(tx["seq"])) if background_under_orfs or outside_id(txs, orfs, name, p)]
        sources += [(outside_id(txs, orfs, name, p) or "background", name, p)
                    for p in rng.choice(free, int(BACKGROUND * len(free)))]
    truth, reads = {}, []
    for comp, name, p in sources:
        seq = txs[name]["seq"]
        length = rng.choice(list(PLANTED))
        ds, ps = zip(*PLANTED[length])
        pos5 = p - rng.choice(ds, p=ps)
        if pos5 < 0 or pos5 + length > len(seq):
            continue
        q = pos5 + 12
        n, s, e = orfs.get(comp, (name, -1, -1))
        own = s <= q < e if comp in orfs else outside_id(txs, orfs, name, q) == comp
        truth[comp if own else "leaked"] = truth.get(comp if own else "leaked", 0) + 1
        fp = seq[pos5:pos5 + length]
        reads.append((f"orf_{len(reads)}_x1", fp, f"{length}M",
                      [(t, i) for t, tx in txs.items() for i in _find_all(tx["seq"], fp)]))
    bam_path = directory / "orfs.bam"
    write_bam(bam_path, txs, reads)

    table = [(i, n, s, e) for i, (n, s, e) in orfs.items() if not i.endswith(":CDS")]
    s1 = orfs["U1_uORF"][1]
    table += [("bad_len", "U1", s1, s1 + 62), ("bad_stop", "U1", s1, s1 + 60),
              ("bad_off", "N1", 300, 330), ("bad_tx", "nope", 0, 30)]
    orf_path = directory / "orfs.tsv"
    orf_path.write_text("ORF_id\tName\tstart\tend\n" + "".join(f"{i}\t{n}\t{s}\t{e}\n" for i, n, s, e in table))
    return dict(fasta=fa, gtf=gtf, bam=bam_path, orfs=orf_path, txs=txs, orf_coords=orfs, truth=truth)
