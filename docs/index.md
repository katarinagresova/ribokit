# ribokit

A Python Ribo-seq toolkit, built milestone by milestone, deterministic from
end to end.

Milestone 1: P-site offset estimation and CDS-level quantification from
transcriptome alignments (`ribokit quant`). Milestone 2: reading-frame
evidence for uORFs and other ORFs besides the CDS (`ribokit orfs` /
`ribokit score`).

## Why

Ribo-seq quantification needs the isoform-length term accounted for, a fit
that actually converges from a fixed, deterministic start, and P-site offsets
estimated per read length from the data rather than assumed. Counting evidence for a
uORF additionally needs a reading frame that isn't thrown away before it's
checked, and a null that accounts for the host CDS's own off-frame noise.
ribokit does all of this, benchmarked against real runs rather than
simulation alone.

## Get started

- **Install and run**: see the [README](https://github.com/katarinagresova/ribokit#readme)
  for setup, the CLI, and output file formats.
- **[Method: quant](method.md)**: how offsets, P-site assignment, and
  quantification work, with interactive examples.
- **[Method: ORF counting](orfs.md)**: ORF types, outside components, the
  frame-weighted EM, and the codon lead score.
- **Using `orfs` and `score`**: a [typical workflow](orfs.md#typical-workflow)
  for two conditions, [how to read the scores](orfs.md#reading-the-results),
  and [what they give on real data](orfs.md#validation-on-real-data).

## Status

Milestone 1 (offsets + CDS quantification) and milestone 2 (ORF/uORF
counting and frame scoring). Later milestones (codon occupancy, the
harringtonine initiation track, ribokit's own ORF caller) are tracked
internally and not part of this release yet.

## License

Apache License 2.0, see [LICENSE](https://github.com/katarinagresova/ribokit/blob/main/LICENSE)
and [NOTICE](https://github.com/katarinagresova/ribokit/blob/main/NOTICE).
