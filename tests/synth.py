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


# Start dataset (make_start_dataset): single-exon transcripts on the + strand, with an elongation and a
# harringtonine library per replicate. Elongation P-sites (expected per nt): each translated ORF's codons at
# its density (reads per codon), a start peak of START_PEAK codons at the annotated starts, a pause, and
# untranslated background per nt (BG_* times the CDS density; BG_TRANSCRIPT without a CDS).
# Harringtonine P-sites: translating ribosomes run off (runoff() of the codon from their ORF's start), the
# background is raised per region (HARR_*), and each planted start adds a kernel of `worth` x density reads
# (KERNEL: share per nt from the start's first nt). Counts are Poisson.
START_PEAK, PAUSE = 5, 15
BG_LEADER, BG_TRAILER, BG_TRANSCRIPT = 0.05, 0.02, 0.05
KERNEL = {0: 0.5, **{3 * c: 0.06 for c in range(1, 6)}, **{3 * c: 0.04 for c in range(6, 11)}}
HARR_LEADER, HARR_TRAILER, HARR_TRANSCRIPT = 4.0, 1.5, 1.5
WORTH = 60     # a planted start: its kernel in codons of its ORF's density
KOZAK = "GCCACC"


def runoff(codon):
    """Harringtonine / elongation of translating ribosomes: 0.5 at the start, 1.4 from codon 300 on."""
    return 0.5 + 0.9 * np.minimum(codon, 300) / 300


def first_stop(seq, s):
    """End of the ORF from s: the first nt of its first in-frame stop, None if there is none."""
    for i in range(s + 3, len(seq) - 2, 3):
        if seq[i:i + 3] in STOPS:
            return i
    return None


def build_start_transcripts(rng, n_fill=30):
    """Transcripts {name: dict(seq, cds_start, cds_end)}, ORFs {id: (name, start, end, density, worth)} (CDSs
    included; worth = the planted start's kernel in codons of density, 0 if not planted) and the pause
    {name: position}. Planted starts sit in a fixed context (KOZAK before, GCC after), so no other start codon
    is within 3 nt of them."""
    rand = lambda n: "".join(rng.choice(list("ACGT"), n))
    codons = lambda n: "".join(rng.choice(SENSE, n))
    stop = lambda: str(rng.choice(STOPS))
    txs, orfs, pauses = {}, {}, {}

    def add(name, leader, n_codons, density, worth=WORTH, trailer=150):
        """A transcript: leader (its last 6 nt KOZAK), then a CDS of n_codons starting ATG GCC, stop, trailer."""
        leader = leader[:-6] + KOZAK
        txs[name] = leader + "ATGGCC" + codons(n_codons - 2) + stop() + rand(trailer)
        orfs[f"{name}:CDS"] = (name, len(leader), len(leader) + 3 * n_codons, density, worth)
        return len(leader)

    for i in range(n_fill):
        add(f"F{i}", rand(int(rng.integers(250, 450))), int(rng.integers(150, 400)), float(rng.uniform(0.5, 3)))
    # S1: a translated ATG uORF and a translated CTG uORF
    lead = rand(54) + KOZAK + "ATGGCC" + codons(11) + stop() + rand(54) + KOZAK
    s = len(lead)
    lead += "CTGGCC" + codons(7) + stop() + rand(60)
    add("S1", lead, 200, 1.5)
    orfs["S1_uATG"] = ("S1", 60, 99, 1.0, WORTH)
    orfs["S1_uCTG"] = ("S1", s, s + 27, 1.0, WORTH)
    # S2: a uoORF, ATG 61 nt before the CDS (another frame), its leader part of sense codons
    lead = rand(84) + KOZAK
    s = len(lead)
    add("S2", lead + "ATGGCC" + codons(18) + rand(7), 200, 1.5)
    orfs["S2_uoORF"] = ("S2", s, first_stop(txs["S2"], s), 1.0, WORTH)
    # S3: an N-terminal extension, CTG 31 codons before the CDS, in frame, no stop
    cs = add("S3", rand(74) + KOZAK + "CTGGCC" + codons(29) + "GCCACC", 200, 1.5)
    orfs["S3_ext"] = ("S3", 80, cs + 600, 1.0, WORTH)
    # S4: a pair, CTG two codons before an ATG uORF start, in frame; the ATG is the start
    lead = rand(74) + KOZAK + "CTGACC"
    s = len(lead)
    lead += "ATGGCC" + codons(9) + stop()
    add("S4", lead + rand(80), 200, 1.5)
    orfs["S4_pairATG"] = ("S4", s, s + 33, 1.0, WORTH)
    orfs["S4_pairCTG"] = ("S4", s - 6, s + 33, 0.0, 0)
    # S5: a weaker start 5 codons after a strong one, in frame +1 (ATG at nt 16 of the first ORF)
    lead = rand(54) + KOZAK
    s = len(lead)
    lead += "ATGGCC" + codons(3) + "CAT" + str(rng.choice([c for c in SENSE if c[0] == "G"])) + codons(13) + stop()
    add("S5", lead + rand(80), 200, 1.5)
    orfs["S5_first"] = ("S5", s, s + 60, 1.0, WORTH)
    orfs["S5_after"] = ("S5", s + 16, first_stop(txs["S5"], s + 16), 0.3, 10)
    # S6: an elongation pause at an in-frame CTG, codon 80 of the CDS: not a start
    cs = add("S6", rand(150), 250, 1.5)
    seq = txs["S6"]
    txs["S6"] = seq[:cs + 240] + "CTG" + seq[cs + 243:]
    pauses["S6"] = cs + 240
    orfs["S6_pause"] = ("S6", cs + 240, cs + 750, 0.0, 0)
    # S7: an untranslated ATG uORF
    add("S7", rand(94) + KOZAK + "ATGGCC" + codons(9) + stop() + rand(100), 200, 1.5)
    orfs["S7_quiet"] = ("S7", 100, 133, 0.0, 0)
    # N1: no CDS, a translated ORF
    txs["N1"] = rand(144) + KOZAK + "ATGGCC" + codons(29) + stop() + rand(200)
    orfs["N1_ORF"] = ("N1", 150, 243, 1.0, WORTH)
    # untranslated candidates: the first ATG or CTG of the first ten fillers' leaders that has a stop
    for i in range(10):
        name = f"F{i}"
        seq, cs = txs[name], orfs[f"{name}:CDS"][1]
        for p in range(cs - 3):
            if seq[p:p + 3] in ("ATG", "CTG") and first_stop(seq, p) is not None:
                orfs[f"{name}_cand"] = (name, p, first_stop(seq, p), 0.0, 0)
                break
    cds = {o[0]: o for i, o in orfs.items() if i.endswith(":CDS")}
    for i, (name, s, e, _, _) in orfs.items():
        assert e is not None and (e - s) % 3 == 0 and txs[name][e:e + 3] in STOPS, i
    return ({name: dict(seq=seq, cds_start=cds[name][1] if name in cds else None,
                        cds_end=cds[name][2] if name in cds else None) for name, seq in txs.items()}, orfs, pauses)


def start_profiles(txs, orfs, pauses):
    """Expected P-sites per nt {name: (elongation, harringtonine)}."""
    out = {}
    for name, tx in txs.items():
        n, cs, ce = len(tx["seq"]), tx["cds_start"], tx["cds_end"]
        e, h = np.zeros(n), np.zeros(n)
        if cs is None:
            e += BG_TRANSCRIPT
            h += BG_TRANSCRIPT * HARR_TRANSCRIPT
        else:
            d = orfs[f"{name}:CDS"][3]
            e[:cs], h[:cs] = BG_LEADER * d, BG_LEADER * d * HARR_LEADER
            e[ce:], h[ce:] = BG_TRAILER * d, BG_TRAILER * d * HARR_TRAILER
            peaks = [(cs, START_PEAK * d)] + ([(pauses[name], PAUSE * d)] if name in pauses else [])
            for p, x in peaks:
                e[p] += x
                h[p] += x * runoff((p - cs) // 3)
        for orf_name, s, end, density, worth in orfs.values():
            if orf_name == name:
                e[s:end:3] += density
                h[s:end:3] += density * runoff(np.arange((end - s) // 3))
                for d, share in KERNEL.items():
                    h[s + d] += worth * density * share
        out[name] = (e, h)
    return out


def make_start_dataset(directory, seed=5, replicates=2):
    """Reference, ORF table (the ORFs but the CDSs) and per replicate an elongation and a harringtonine BAM
    (libs[(kind, rep)])."""
    rng = np.random.default_rng(seed)
    txs, orfs, pauses = build_start_transcripts(rng)
    genome, pos = [], 0
    for name, tx in txs.items():
        genome += ["".join(rng.choice(list("ACGT"), 150)), tx["seq"]]
        tx.update(gene=f"g_{name}", strand="+", gexons=[(pos + 150, pos + 150 + len(tx["seq"]))])
        pos += 150 + len(tx["seq"])
    fa, gtf = write_reference(directory, "".join(genome) + "A" * 150, txs)
    prof = start_profiles(txs, orfs, pauses)
    libs = {}
    for rep in range(1, replicates + 1):
        for k, kind in enumerate(("elongation", "harringtonine")):
            reads = []
            for name, tx in txs.items():
                seq = tx["seq"]
                for p in np.repeat(np.arange(len(seq)), rng.poisson(prof[name][k])):
                    length = int(rng.choice(list(PLANTED)))
                    ds, ps = zip(*PLANTED[length])
                    pos5 = p - int(rng.choice(ds, p=ps))
                    if 0 <= pos5 <= len(seq) - length:
                        reads.append((f"{kind}{rep}_{len(reads)}_x1", seq[pos5:pos5 + length], f"{length}M",
                                      [(name, pos5)]))
            libs[kind, rep] = directory / f"{kind}{rep}.bam"
            write_bam(libs[kind, rep], txs, reads)
    table = [(i, n, s, e) for i, (n, s, e, _, _) in orfs.items() if not i.endswith(":CDS")]
    orf_path = directory / "starts_orfs.tsv"
    orf_path.write_text("ORF_id\tName\tstart\tend\n" + "".join(f"{i}\t{n}\t{s}\t{e}\n" for i, n, s, e in table))
    return dict(fasta=fa, gtf=gtf, orfs=orf_path, libs=libs, txs=txs, orf_coords=orfs, pauses=pauses)
