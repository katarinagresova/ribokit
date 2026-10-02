# ribokit

Ribo-seq analysis in Python. Milestone 1 (this version) does one job: P-site
offset estimation and CDS quantification from transcriptome alignments,
deterministic from end to end.

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
| `<prefix>.quant.tsv` | `Name Length EffectiveLength ritpm NumReads`, one row per CDS, sorted by Name. `Length` = CDS without the stop codon. `NumReads` = expected reads; they sum to the assigned reads. `ritpm` = reads per nt, scaled to sum to 1e6. `NA` where no read is compatible. |
| `<prefix>.offsets.tsv` | per (length, phase): `offset` (nt from the 5' end to the P-site codon), `support` (reads spanning a start or stop), `z` (how firmly the data pin this phase's offset: score lead of the chosen window over the best window that gives this phase another offset, divided by the square root of the reads the two disagree on; below about 4 = weakly determined), `reads` |
| `<prefix>.ties.tsv` | CDSs that no read tells apart (e.g. identical paralogs); their split comes from the EM's even start, not from the data |
| `<prefix>.stats.tsv` | reads left after each filter, EM steps and log-likelihood |
| `<prefix>.psites.tsv` | `read Name psite length`, one row per alignment assigned to a CDS: `psite` is the 0-based transcript position of the P-site's first nt |

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

## License

Apache License 2.0, see [LICENSE](LICENSE) and [NOTICE](NOTICE).
