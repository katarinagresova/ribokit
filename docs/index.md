# ribokit

A Python Ribo-seq toolkit, built milestone by milestone.

Milestone 1 does one job: P-site offset estimation and CDS-level
quantification from transcriptome alignments, deterministic from end to
end.

## Why

Ribo-seq quantification needs the isoform-length term accounted for, a fit
that actually converges from a fixed, deterministic start, and read-length
windows chosen from the data rather than assumed. ribokit does all three,
benchmarked against real runs rather than simulation alone.

## Get started

- **Install and run**: see the [README](https://github.com/katarinagresova/ribokit#readme)
  for setup, the CLI, and output file formats.
- **[Method](method.md)**: how offsets, P-site assignment, and quantification
  work, with interactive examples.

## Status

Milestone 1 only: offsets + CDS quantification. Later milestones (shared
core, elongation track, harringtonine initiation track) are tracked
internally and not part of this release yet.

## License

Apache License 2.0, see [LICENSE](https://github.com/katarinagresova/ribokit/blob/main/LICENSE)
and [NOTICE](https://github.com/katarinagresova/ribokit/blob/main/NOTICE).
