import filecmp

import pandas as pd
import pysam
import pytest

from ribokit import orfs, quant
from ribokit.bam import _footprint_length

from synth import PLANTED, make_dataset

WINDOW = (28, 30)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("synth"))


@pytest.fixture(scope="module")
def run(data, tmp_path_factory):
    prefix = tmp_path_factory.mktemp("orfs") / "s"
    frames, off, stats = orfs.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(prefix))
    return dict(prefix=prefix, frames=frames.set_index("length"), offsets=off, stats=stats)


def test_offsets_as_quant(run, data, tmp_path):
    quant.quantify(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(tmp_path / "q"))
    assert filecmp.cmp(f"{run['prefix']}.offsets.tsv", tmp_path / "q.offsets.tsv", shallow=False)


def test_psites_one_offset_per_length(run, data):
    # every forward alignment in the window, on any reference (UTRs, gN without a CDS, the
    # reference missing from the GTF), P-site = 5' end + the planted phase-0 offset (12)
    expected = []
    with pysam.AlignmentFile(str(data["bam"])) as f:
        for a in f:
            fp = _footprint_length(a.cigartuples)
            if fp is not None and WINDOW[0] <= fp <= WINDOW[1]:
                expected.append((a.query_name, a.reference_name, a.reference_start + 12, fp))
    got = pd.read_csv(f"{run['prefix']}.psites.tsv", sep="\t")
    assert list(got.itertuples(index=False, name=None)) == expected
    assert run["stats"]["reads_with_psite"] == run["stats"]["reads_in_length_window"]


def test_frames_match_planted(run):
    # a read whose true offset is d has its P-site at true + 12 - d: frame (12 - d) mod 3
    for length, ds in PLANTED.items():
        row = run["frames"].loc[length]
        assert row["offset"] == 12
        expected = [sum(p for d, p in ds if (12 - d) % 3 == f) for f in range(3)]
        assert [row["frame0"], row["frame1"], row["frame2"]] == pytest.approx(expected, abs=0.03)
        assert row["reads"] > 1000


def test_deterministic_and_offsets_reusable(run, data, tmp_path):
    again = tmp_path / "again"
    orfs.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(again))
    reused = tmp_path / "reused"
    orfs.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(reused),
             offsets_path=f"{run['prefix']}.offsets.tsv")
    for suffix in ("frames.tsv", "offsets.tsv", "psites.tsv", "stats.tsv"):
        assert filecmp.cmp(f"{run['prefix']}.{suffix}", f"{again}.{suffix}", shallow=False), suffix
    for suffix in ("frames.tsv", "psites.tsv"):
        assert filecmp.cmp(f"{run['prefix']}.{suffix}", f"{reused}.{suffix}", shallow=False), suffix
