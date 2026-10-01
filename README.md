# ribokit

Ribo-seq analysis in Python. Milestone 1 (this version) does one job: P-site
offsets and CDS quantification from transcriptome alignments, as a drop-in
replacement for the RiboStan steps in wf-eIF-deltaTE. See PLAN.md for the
review of the design record, the scope, and what comes later.

## Install

```bash
mamba env create -p .env -f environment.yml
.env/bin/pip install --no-deps -e .
.env/bin/python -m pytest          # 9 tests on synthetic data with known truth
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
| `<prefix>.quant.tsv` | `Name Length EffectiveLength ritpm NumReads`, one row per CDS, sorted by Name. `Length` = CDS without the stop codon. `NumReads` = expected reads; they sum to the assigned reads. `ritpm` = reads per nt, scaled to sum to 1e6. `NA` where no read is compatible. Same columns as RiboStan's `*_morf_quant.tsv`. |
| `<prefix>.offsets.tsv` | per (length, phase): `offset` (nt from the 5' end to the P-site codon), `support` (reads spanning a start or stop), `margin` (score lead of the chosen window over the next, per read of support; near 0 = weakly determined), `reads` |
| `<prefix>.ties.tsv` | CDSs that no read tells apart (e.g. identical paralogs); their split comes from the EM's even start, not from the data |
| `<prefix>.stats.tsv` | reads left after each filter, EM iterations and log-likelihood |

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
   P(read | CDS) = 1/length. It starts deterministically and runs until no
   expected count moves by more than 1e-3 reads. It fails if it does not get
   there.

## License

Apache License 2.0, see [LICENSE](LICENSE). ribokit is a reimplementation, not
a port: it contains no RiboStan code (RiboStan is GPL-3).
