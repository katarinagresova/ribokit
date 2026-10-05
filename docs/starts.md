# Start sites from harringtonine

`ribokit starts` finds translation start sites from harringtonine (or
lactimidomycin) Ribo-seq libraries. It compares each harringtonine library
with an elongation library of the same condition and replicate, tests given
starts or every start codon of the transcriptome, calls starts, and writes
their ORFs and a catalogue of uORFs and uoORFs.

## Background

**Harringtonine stops ribosomes at the start codon.** Harringtonine binds the
large ribosomal subunit and blocks the first steps of elongation, but not
initiation ([Fresno et al. 1977](#references)). Ribosomes that start after the
drug is added stay on the start codon; ribosomes that were already elongating
run off the mRNA. After a few minutes, footprints pile up on start codons
([Ingolia et al. 2011](#references)). Lactimidomycin acts in the same way
([Lee et al. 2012](#references)).

**A peak in a harringtonine library is not yet a start.** Three things also
make harringtonine footprints:

- **Run-off is not complete.** Ribosomes far into a long CDS have not run off
  yet, so the harringtonine / elongation ratio grows along the CDS. Leaders
  have a higher ratio than CDSs.
- **Codons 1-10 after a start are raised.** Ribosomes that made a few peptide
  bonds before the drug acted stay there: 5-8x at codons 1-5, about 3x at
  codons 6-10, against the CDS.
- **Pauses.** A codon where elongating ribosomes stall has a peak in every
  library, with or without the drug.

Thus `starts` does not look for peaks in the harringtonine library alone. It
asks if a start codon has more harringtonine P-sites than the elongation
P-sites at the same position predict.

## Inputs

`starts` reads the outputs of `ribokit orfs` runs
([ORF counting](orfs.md)), one per library, and a library sheet:

```text
pool        harringtonine          elongation
minusAux    out/harr_minus_rep1    out/elong_minus_rep1
minusAux    out/harr_minus_rep2    out/elong_minus_rep2
plusAux     out/harr_plus_rep1     out/elong_plus_rep1
plusAux     out/harr_plus_rep2     out/elong_plus_rep2
```

Each row is one replicate: a harringtonine library and its matched elongation
library, given as out prefixes of `orfs`. A **pool** is the replicates of one
condition. Several pools can be given; each keeps its own background, and
their evidence adds up.

The libraries of a pool must have the same P-site offsets, or their P-sites
are not comparable. Give each `orfs` run of a pool the offsets of
`ribokit quant` on the pool's harringtonine BAMs, merged into one
(`--offsets`). The libraries also need the same transcripts: give them the
same GTF. Transcripts that are not in every `orfs` run are left out.

```bash
samtools cat -o harr_minus.bam harr_minus_rep1.bam harr_minus_rep2.bam
ribokit quant --bam harr_minus.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/minus
for lib in harr_minus_rep1 harr_minus_rep2 elong_minus_rep1 elong_minus_rep2; do
    ribokit orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-30 --offsets out/minus.offsets.tsv --out-prefix out/$lib
done
# the same for the plusAux pool, with its own offsets
ribokit starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
```

## Method

### 1. Start lengths and the kernel

For each read length, `starts` counts the harringtonine P-sites around the
annotated start codons. The **start lengths** are the longest run of
consecutive lengths whose P-sites peak on the first nt of the start codon in
every harringtonine library (`--start-lengths` sets them). Only their P-sites
count, one per alignment, in both kinds of library. Shorter and longer reads
often peak 1-3 nt away, because their offsets do not fit the start.

The **kernel** is the shape of a start: per length, the harringtonine P-sites
minus the elongation P-sites (scaled on codons 11-20), as a share of nt -15
to +32. `kernel.tsv` has it per pool.

### 2. The window

A start's evidence is the P-sites in its **window**: nt -1 to +2 of its codon
(the codon and the nt before it). The kernel puts about half of a start's
P-sites on nt 0, a few % on nt -1, and the rest on codons 1-10. A wider
window, such as the codon +-1 codon, holds a few more of them, but it also
holds nt 0 of the codon before or after. Then two start codons one codon
apart get the same reads, and the scan cannot tell which one is used. With
the narrow window, the start codon wins against its neighbours in 97-98% of
the annotated starts.

### 3. The background

Without a start, the harringtonine P-sites of a window follow the elongation
P-sites there, times a factor that depends on where the window is:

$$
\text{expected } h \;=\; s_t \, f_b \, \lambda .
$$

- $\lambda$ is the window's elongation rate (step 4).
- $f_b$ is a factor per **bin** of region and distance: leaders by nt
  before the CDS start (0, 30, 100, 300 and more), the CDS by codon from its
  start (0, 20, 50, 100, 200, 400, 700, 1,000 and more), trailers by nt after
  the stop codon, and transcripts without an annotated CDS.
- $s_t$ is a scale per transcript: harringtonine against elongation reads
  differ between transcripts.

The factors are fit on **null codons**, by a Poisson fit with alternating
updates. A null codon is a position with no start codon (ATG or one of its nine
near-cognates) within 8 nt, outside the window and codons 1-10 of an annotated
start, and with mostly unique reads. (Processed pseudogenes share their
parent's reads, start peak included.) The scale of a transcript comes from all
its positions but the windows and codons 1-10 of annotated starts.
`factors.tsv` has the factors.

### 4. The test

Write $R = s_t f_b$, $h$ and $e$ for the harringtonine and elongation P-sites
in the window. The window's elongation rate has a gamma prior. Its mean
$\lambda_0$ is the elongation density of the window's nt in their transcript:
per zone (leader, CDS frame 0, 1 and 2, trailer), as a window of 4 nt holds
one or two frame-0 nt of a CDS. Its shape $a$ is fit per region on the null
windows. Then

$$
h \sim \text{NB}\Big(\text{mean } m = R \max\!\big(e,\; \tfrac{a + e}{a/\lambda_0 + 1}\big),\;
\text{shape } a + e\Big),
$$

with an extra overdispersion $\phi$ (variance $+\,\phi m^2$). The posterior
mean lifts a window whose elongation reads are few by chance; the maximum keeps
an elongation peak (a pause) from being shrunk, so a pause is not called.

- **Pools** add their means and variances.
- **$\phi$** is calibrated per region and bin of expected reads on the null
  windows (one per 4 nt): the smallest $\phi$ that puts at most 1% of them at
  $p \le 0.01$.
- $p = P(H \ge h)$; `enrichment` is $h / m$.
- $q$ is Benjamini-Hochberg per **class**: the annotated starts, and the other
  starts per region (leader, CDS, trailer, transcript without a CDS). With
  `--scan`, BH counts every start codon of the class: the ones without
  harringtonine reads are tests with $p = 1$.
- `typical_p` is the $p$ the window would have with the reads of a typical
  used start there: the median enrichment of the annotated starts with 3 or
  more expected reads. If it is large, the start cannot be called.
- `starts_stats.tsv` gives, per class, the rate of null windows at the $p$ of
  the last call, and the false calls that rate gives over the class's tests.

### 5. The scan

With `--scan`, the candidates are every ATG, CTG, GTG, TTG, AAG, ACG, AGG,
ATA, ATC and ATT with harringtonine reads in its window, and the annotated
starts. Starts near each other share reads, so they are called in order:

1. The codons 1-10 of an annotated start are not called: harringtonine raises
   them when the start is used, and the test cannot tell them from a start.
2. The candidates are taken by $p$, smallest first. A candidate with
   $q < 0.05$ is **called** unless a called start stops it.
3. A called start stops the candidates within 8 nt of it and in its codons
   1-10. `stopped_by` names the start that stops a candidate.

Thus of two starts one or two codons apart, only the stronger is called.

### 6. ORFs and the catalogue

Each called start runs to its first in-frame stop codon. Its type comes from
the coordinates against the annotated CDS: `uORF` (ends at or before the CDS
start), `uoORF` (starts before the CDS, out of frame, ends in it), `extension`
(starts before the CDS, in frame), `CDS`, `truncation` (inside the CDS, in
frame), `internal` (inside the CDS, out of frame), `dORF` (after the CDS) or
`other` (no annotated CDS). `start_orfs.tsv` has them, but the annotated CDSs.
`catalogue.tsv` has the uORFs and uoORFs, and every annotated CDS with its
start evidence. Both are ORF tables for `ribokit orfs --orfs`: count the
catalogue's ORFs in other libraries, and score their frames with
`ribokit score`.

## Reading the results

- **A call is start evidence, not translation evidence.** It says that
  ribokit found ribosomes stopped on that codon by harringtonine. Whether
  ribosomes go on to elongate is a frame question: `orfs` and `score` on the
  ORF table.
- **Power depends on the elongation reads.** A start with few elongation reads
  in its window has a wide posterior for its elongation rate, so it needs many
  harringtonine reads. Look at `expected` and `typical_p`.
- **Leaders need more.** Leaders have a higher harringtonine / elongation
  ratio than CDSs near their start, so a leader start needs more
  harringtonine reads than a CDS start with the same elongation reads.
- **Trailers.** Harringtonine libraries can hold RNA fragments that the
  elongation libraries do not have, mostly in trailers and of short reads.
  `starts_stats.tsv` gives many expected false calls in trailers; use those
  calls with care.
- **Multimappers.** Each alignment counts once, so paralogs and pseudogenes
  share a start peak. `multimapped` gives the share of the window's reads
  with another alignment.

## Validation on real data

HCT116, 4 h, two pools of three replicates (one depletion line, without and
with auxin): six harringtonine libraries and their six elongation libraries.

- **Start lengths**: 23-29 nt. The kernel holds 44-51% of its P-sites at nt 0
  and 3-4% at nt -1.
- **Background**: leader factors are 2.4-2.9 near the CDS start and 1.0-1.4
  300 nt or more before it; CDS factors are 0.16-0.22 at codons 20-200 and
  reach 1 at codon 1,000 (all against codons 1,000 and more).
- **Null windows** (held-out transcripts): 0.17-0.19x alpha at $p \le 0.01$
  and 0.9-1.2x alpha at $p \le 0.001$ in leaders.
- **Annotated starts** (scan): 60.5% called; 29% of those with less than 0.3
  expected reads in the window, 82% with 1-3, 93-95% with 3 or more. Near an
  annotated start, the best call is on the annotated codon in 99.4%.
- **Swaps**: with an elongation library as the "harringtonine" library of
  another replicate (replicate 3 against 2, in each pool), 6 of 12,294
  annotated starts are called. With replicate 1, whose short footprints fit
  other offsets than the pool's, 410-417 are: give each pool offsets that fit
  all its libraries.
- **Start context**: called leader starts have a purine at -3 in 58-75% and a
  G at +4 in 40-60%; start codons in leaders that are not called, at the
  same depth, 50-58% and 22-37% ([Kozak 1986](#references)). The sequence
  plays no part in the test.

## References

- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289-300.
- Fresno M, Jiménez A, Vázquez D (1977). Inhibition of translation in
  eukaryotic systems by harringtonine. *Eur J Biochem* 72:323-330.
- Ingolia NT, Lareau LF, Weissman JS (2011). Ribosome profiling of mouse
  embryonic stem cells reveals the complexity and dynamics of mammalian
  proteomes. *Cell* 147:789-802.
- Kozak M (1986). Point mutations define a sequence flanking the AUG initiator
  codon that modulates translation by eukaryotic ribosomes. *Cell* 44:283-292.
- Lee S, Liu B, Lee S, Huang SX, Shen B, Qian SB (2012). Global mapping of
  translation initiation sites in mammalian cells at single-nucleotide
  resolution. *PNAS* 109:E2424-E2432.
