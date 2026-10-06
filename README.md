# ribokit

Ribo-seq analysis in Python, deterministic from end to end: P-site offset
estimation and CDS quantification from transcriptome alignments
(`ribokit quant`), reading-frame evidence for uORFs and other non-CDS ORFs
(`ribokit orfs` / `ribokit score`), and translation start sites from
harringtonine libraries (`ribokit starts`).

## Install

```bash
mamba env create -p .env -f environment.yml
.env/bin/pip install --no-deps -e .
.env/bin/python -m pytest          # 43 tests on synthetic data with known truth
```

## Use

```bash
ribokit quant --bam sample.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/sample.human
# spike-in, reusing the human offsets (or leave --offsets out to estimate them from the spike-in reads)
ribokit quant --bam sample.spikein.bam --gtf yeast.gtf --fasta yeast.fa \
    --read-lengths 18-30 --offsets out/sample.human.offsets.tsv --out-prefix out/sample.yeast
```

Inputs:
- `--bam`: alignments to transcripts; reference names are transcript ids. Text
  after a `|` is ignored.
- `--gtf` / `--fasta`: the annotation and its genome. Transcript lengths must
  match the BAM header.
- `--read-lengths`: footprint lengths to use (required).

| Output | Content |
|---|---|
| `<prefix>.quant.tsv` | `Name Length EffectiveLength ritpm NumReads`, one row per CDS, sorted by Name. `Length` = CDS without the stop codon. `NumReads` = expected reads; they sum to the assigned reads. `ritpm` = reads per nt, scaled to sum to 1e6. `NA` where no read is compatible. |
| `<prefix>.offsets.tsv` | per (length, phase): `offset` (nt from the 5' end to the P-site codon), `support` (reads covering the first nt of the start codon or of the last sense codon; a read with n alignments counts 1/n), `z` (how firmly the data pin this phase's offset: score lead of the chosen window over the best window that gives this phase another offset, divided by the square root of the reads the two disagree on; below about 4 = weakly determined), `reads` |
| `<prefix>.ties.tsv` | `Name tie_group`: CDSs that no read tells apart (e.g. identical paralogs). Their split comes from the model, not from the data: equally long CDSs keep the EM's even start, otherwise the shortest gets almost all the reads. |
| `<prefix>.stats.tsv` | reads left after each filter, EM steps and log-likelihood |
| `<prefix>.psites.tsv` | `read Name psite length`, one row per alignment assigned to a CDS: `psite` is the 0-based transcript position of the P-site's first nt |

### ORFs other than the CDS: `orfs` and `score`

`ribokit orfs` counts reads on a table of ORFs — upstream ORFs (uORFs),
ORFs overlapping the CDS out of frame (uoORFs), and others — alongside the
annotated CDS, and keeps every read's reading frame instead of snapping it to
the CDS. `ribokit score` turns those per-codon counts into frame evidence for
each ORF, pooled over one or more `orfs` runs (e.g. replicate libraries):

```bash
ribokit orfs --bam sample_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --orfs candidates.tsv --out-prefix out/sample_rep1
# the same for sample_rep2, with --offsets out/sample_rep1.offsets.tsv
ribokit score --orfs-prefix out/sample_rep1 out/sample_rep2 --out-prefix out/pooled
```

`--orfs` is a TSV `ORF_id Name start end` (0-based transcript coordinates,
`start` = first nt of the start codon, `end` = one past the last sense codon,
same convention as the CDS). The annotated CDSs are added automatically.

To compare conditions, give every library the same GTF, ORF table and
offsets table, and score all libraries together:
[Typical workflow](docs/orfs.md#typical-workflow). What the scores do and do
not say: [Reading the results](docs/orfs.md#reading-the-results).

| Output | Content |
|---|---|
| `<prefix>.orfs.tsv` | `ORF_id Name type start end start_codon Length NumReads`, one row per ORF (from `--orfs`, plus the annotated CDSs) and per outside component (`<tx>:leader`, `<tx>:trailer`, `<tx>:transcript`: the positions in no ORF). `type` is `CDS`, `uORF` (ends at or before the CDS start), `uoORF` (starts upstream of the CDS, out of frame with it) or `other` for an ORF, and `leader`, `trailer` or `transcript` for an outside component. |
| `<prefix>.offsets.tsv` | the same table as `quant`'s (the same BAM and read lengths give the same offsets); either command's `--offsets` accepts it. `orfs` uses only each length's phase-0 offset. |
| `<prefix>.frames.tsv` | per read length: `offset` (the phase-0 offset, not a window), `reads`, and `frame0 frame1 frame2` — the share of P-sites in each frame of the CDS, from CDS interiors. |
| `<prefix>.codons.tsv` | `ORF_id codon length frame0 frame1 frame2`: P-sites per codon of every ORF but the annotated CDSs, by length and frame — `score`'s input. |
| `<prefix>.ties.tsv`, `<prefix>.stats.tsv`, `<prefix>.psites.tsv` | as for `quant`, but over every ORF and outside component: `ties.tsv` is `ORF_id tie_group`, and `psites.tsv` keeps every P-site, not only those assigned to a CDS. |
| `<prefix>.scores.tsv` | one row per scored ORF: `ORF_id Name type start end codons codons_with_reads reads in_frame_share leads expected z p min_p q` — `p`/`q` are the frame-evidence significance and its BH correction; `min_p` is the best `p` the ORF's length and reads could reach. |
| `<prefix>.decoys.tsv` | each ORF shifted by +1 and +2 nt and scored the same way, for null calibration. |
| `<prefix>.score_stats.tsv` | the read lengths used and how many ORFs/decoys were scored. |

### Start sites: `starts`

`ribokit starts` finds translation start sites in harringtonine (or
lactimidomycin) libraries. Each harringtonine library is compared with a
matched elongation library (same condition and replicate), so that a pause,
which both have, is not taken for a start. It reads the outputs of
`ribokit orfs` runs:

```bash
# one offset table per pool: quant on the merged harringtonine BAMs
ribokit quant --bam harr_merged.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/pool
# orfs on every library, harringtonine and elongation, with that table
ribokit orfs --bam harr_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --offsets out/pool.offsets.tsv --out-prefix out/harr_rep1
# ... harr_rep2, elong_rep1, elong_rep2 the same way
ribokit starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
```

`--libraries` is a TSV `pool harringtonine elongation`, one row per
replicate, each an out prefix of an `orfs` run. A pool is a set of
replicates of one condition. Several pools can be given, and their evidence
adds up. `--scan` tests every ATG and near-cognate start codon (CTG, GTG,
TTG, AAG, ACG, AGG, ATA, ATC, ATT) with harringtonine reads, and calls the
starts with BH q < 0.05; `--start-codons ATG,CTG` tests only these codons and
`--fdr 0.01` sets the cut. `--orfs table.tsv` tests the starts of an ORF
table instead (columns `ORF_id Name start`). The annotated starts are always
tested. How it works and how to read it: [docs/starts.md](docs/starts.md).

| Output | Content |
|---|---|
| `<prefix>.starts.tsv` | one row per tested start: `Name start codon region frame ORF_id harringtonine elongation expected enrichment p typical_p q multimapped`. `harringtonine` and `elongation` are the P-sites in the start's window (nt -1 to +2 of the start codon), summed over the replicates and pools; `expected` is the harringtonine P-sites that the elongation profile predicts there; `q` is BH per class (annotated starts; other starts per region: leader, CDS, trailer, transcript without a CDS); `typical_p` is the p the window would have with the reads of a typical used start; `multimapped` is the share of the window's reads with another alignment. With `--scan` also `called stopped_by end type`. |
| `<prefix>.start_orfs.tsv` | with `--scan`: the ORF of each called start, to its first in-frame stop, but the annotated CDSs: `ORF_id Name start end type start_codon` and the start evidence. `ORF_id` is `<transcript>_<1-based start>`. `type` is `uORF`, `uoORF`, `extension` (in frame upstream of the CDS), `truncation` (in frame inside it), `internal` (out of frame inside it), `dORF` (after it) or `other` (no annotated CDS). An ORF table for `ribokit orfs --orfs`. |
| `<prefix>.catalogue.tsv` | with `--scan`: the uORFs and uoORFs of the called starts, and every annotated CDS with its start evidence; the columns of `start_orfs.tsv`. An ORF table for `ribokit orfs --orfs`. |
| `<prefix>.kernel.tsv` | per pool, read length and nt from the annotated start (-15 to +32): harringtonine and elongation P-sites, and the kernel (harringtonine minus scaled elongation, as a share). |
| `<prefix>.factors.tsv` | per pool: the harringtonine / elongation factor of each bin of region and distance. |
| `<prefix>.starts_stats.tsv` | the start lengths, P-site totals, the prior and overdispersion per region, and per class the tests, the calls and the false calls expected from null windows. |

## Method, in short

1. **Annotation.** One CDS per transcript, in transcript coordinates, without
   the stop codon. It is dropped if its length is not a multiple of 3 or no
   stop follows it. Non-ATG starts are kept.
2. **Reads.** Mapped, forward-strand, primary and secondary alignments. A 5'
   soft clip counts as untemplated nucleotides: the 5' end is the first
   aligned base and the length excludes the clip. Alignments with
   indels/splices are dropped. Each read name counts once.
3. **Offsets.** CDS inclusion: for each read length, the three consecutive
   offsets (one per phase) that put the most reads' P-sites inside CDSs.
4. **Assignment.** A read is compatible with a CDS if its P-site falls inside
   it.
5. **EM.** Expectation-maximisation over equivalence classes, with
   P(read | CDS) = 1/length, accelerated with SQUAREM. It starts
   deterministically and runs until an EM step moves no expected count by
   more than 1e-3 reads. It fails if it does not get there.

`orfs` and `score` add:

6. **ORF table.** uORFs, uoORFs and other ORFs besides the CDS, in transcript
   coordinates; dropped if not a multiple of 3 or no stop follows. Positions
   in no ORF form per-transcript leader/trailer/whole-transcript components.
   Alignments to BAM references that are not in the GTF are dropped.
7. **Frame-true P-sites.** One offset per read length (no phase snapping), so
   every P-site keeps its reading frame, plus a per-length frame profile from
   CDS interiors.
8. **EM with a frame term.** As above, but P(read | ORF) also depends on how
   often the read length lands in that ORF's frame, so reads shared with the
   CDS split by frame as well as by density.
9. **Codon lead score.** Per ORF, each codon (its start excluded) votes
   whether more of its reads are in the ORF's frame than in either other
   frame; the score is how unlikely that many leads are if the ORF were not
   translated, against a null that accounts for a host CDS's own frame bleed.

`starts` adds:

10. **Start lengths.** The read lengths whose harringtonine P-sites peak on
    the first nt of the annotated start codons in every library.
11. **Background.** The elongation P-sites times a harringtonine factor per
    region (leader, CDS, trailer) and distance from the CDS start, and a
    scale per transcript, fit on null codons: positions with no start codon
    within 8 nt.
12. **The test.** Per start, the harringtonine P-sites in nt -1 to +2 of its
    codon against the expected, negative binomial, with a gamma prior on the
    window's elongation rate and an overdispersion calibrated on the null
    codons.
13. **The scan.** Start codons with the smallest p first; a start is called
    if it passes BH and no called start within 8 nt or 10 codons upstream
    stops it. Each called start gives an ORF to its first in-frame stop.

Full method, with worked examples: [docs/method.md](docs/method.md) (`quant`),
[docs/orfs.md](docs/orfs.md) (`orfs` / `score`) and
[docs/starts.md](docs/starts.md) (`starts`).

## License

Apache License 2.0, see [LICENSE](LICENSE) and [NOTICE](NOTICE).
