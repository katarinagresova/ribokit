import filecmp

import numpy as np
import pandas as pd
import pysam
import pytest

from ribokit import orfs, quant
from ribokit.bam import _footprint_length

from synth import PLANTED, make_dataset, make_orf_dataset

WINDOW = (28, 30)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("synth"))


@pytest.fixture(scope="module")
def run(data, tmp_path_factory):
    prefix = tmp_path_factory.mktemp("orfs") / "s"
    table, frames, off, stats = orfs.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(prefix))
    return dict(prefix=prefix, table=table, frames=frames.set_index("length"), offsets=off, stats=stats)


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
    for suffix in ("orfs.tsv", "frames.tsv", "offsets.tsv", "ties.tsv", "psites.tsv", "stats.tsv"):
        assert filecmp.cmp(f"{run['prefix']}.{suffix}", f"{again}.{suffix}", shallow=False), suffix
    for suffix in ("orfs.tsv", "frames.tsv", "psites.tsv"):
        assert filecmp.cmp(f"{run['prefix']}.{suffix}", f"{reused}.{suffix}", shallow=False), suffix


@pytest.fixture(scope="module")
def orf_data(tmp_path_factory):
    return make_orf_dataset(tmp_path_factory.mktemp("orfsynth"))


@pytest.fixture(scope="module")
def orf_run(orf_data, tmp_path_factory):
    prefix = tmp_path_factory.mktemp("orfrun") / "s"
    table, frames, off, stats = orfs.run(str(orf_data["bam"]), str(orf_data["gtf"]), str(orf_data["fasta"]), WINDOW,
                                         str(prefix), orfs_path=str(orf_data["orfs"]))
    return dict(prefix=prefix, table=table.set_index("ORF_id"), stats=stats)


def test_orf_table(orf_run, orf_data):
    t, s = orf_run["table"], orf_run["stats"]
    assert (s["orfs_dropped_transcript_not_in_gtf"], s["orfs_dropped_not_multiple_of_3"],
            s["orfs_dropped_off_transcript"], s["orfs_dropped_no_stop"]) == (1, 1, 1, 1)
    assert s["orfs_cds_added"] == 3            # U1, U3 and U4; U2's is listed as U2_main
    orf_types = {"U1_uORF": "uORF", "U1_cand": "uORF", "U1_CUG": "uORF", "U2_uoORF": "uoORF", "U3_uoORF": "uoORF",
                 "U4_uoORF": "uoORF", "U1:CDS": "CDS", "U2_main": "CDS", "U3:CDS": "CDS", "U4:CDS": "CDS",
                 "N1_ORF": "other"}
    assert t.loc[t["type"].isin(orfs.ORF_TYPES), "type"].to_dict() == orf_types
    assert t.at["U1_CUG", "start_codon"] == "CTG"
    for i, (name, start, end) in orf_data["orf_coords"].items():
        assert tuple(t.loc[i, ["Name", "start", "end"]]) == (name, start, end)
    outside = {"U1:leader", "U1:trailer", "U2:leader", "U2:trailer", "U3:leader", "U3:trailer", "U4:leader",
               "U4:trailer", "P1:transcript", "N1:transcript", "unannotated:transcript"}
    assert set(t.index[~t["type"].isin(orfs.ORF_TYPES)]) == outside
    # a leader's positions are those in no ORF: U1's 277 nt hold uORFs of 63, 48 and 33 nt
    assert t.at["U1:leader", "Length"] == 277 - 63 - 48 - 33


def test_orf_counts_match_truth(orf_run, orf_data):
    n, truth = orf_run["table"]["NumReads"], orf_data["truth"]
    for i in ("U1_uORF", "U1_CUG", "N1_ORF", "U1:CDS"):
        assert abs(n[i] - truth[i]) / truth[i] < 0.02, i
    assert pd.isna(n["U1_cand"])                 # untranslated, no background under ORFs
    # uoORF and CDS split their overlap by frame: shares as in the cassette test
    for uo, cds in (("U2_uoORF", "U2_main"), ("U3_uoORF", "U3:CDS")):
        assert abs(n[uo] / (n[uo] + n[cds]) - truth[uo] / (truth[uo] + truth[cds])) < 0.02, uo
        assert abs(n[uo] + n[cds] - truth[uo] - truth[cds]) / (truth[uo] + truth[cds]) < 0.02, uo
    # P1 holds part of U1's CDS: shared reads split between U1's CDS and P1's outside component
    assert abs(n["P1:transcript"] + n["U1:CDS"] - truth["P1:transcript"] - truth["U1:CDS"]) < 0.01 * truth["U1:CDS"]
    assert n.sum() == pytest.approx(orf_run["stats"]["reads_with_psite"])
    assert orf_run["stats"]["tie_groups"] == 0


def test_orf_run_deterministic(orf_run, orf_data, tmp_path):
    orfs.run(str(orf_data["bam"]), str(orf_data["gtf"]), str(orf_data["fasta"]), WINDOW, str(tmp_path / "again"),
             orfs_path=str(orf_data["orfs"]))
    for suffix in ("orfs.tsv", "ties.tsv", "stats.tsv"):
        assert filecmp.cmp(f"{orf_run['prefix']}.{suffix}", tmp_path / f"again.{suffix}", shallow=False), suffix


def test_frame_term_returns_overlap_reads_to_cds(orf_run, orf_data, tmp_path, monkeypatch):
    # U4_uoORF is untranslated, but its leader part holds the CDS start peak's reads that are 1-2 nt
    # long at the 5' end, so the EM takes it for translated. Split by density alone, the overlap
    # gives it ~90 reads; the frame term returns >= 15 of them to the CDS (21-35 over 10 seeds).
    # It keeps ~60: CDS reads in frame 1 are in its frame 0.
    monkeypatch.setattr(orfs, "frame_weights", lambda frames: np.ones((len(frames), 3)))
    t, *_ = orfs.run(str(orf_data["bam"]), str(orf_data["gtf"]), str(orf_data["fasta"]), WINDOW,
                     str(tmp_path / "noframe"), orfs_path=str(orf_data["orfs"]))
    density_only = t.set_index("ORF_id")["NumReads"]
    n = orf_run["table"]["NumReads"]
    assert n["U4:CDS"] - density_only["U4:CDS"] >= 15
    assert abs(n["U4:CDS"] - orf_data["truth"]["U4:CDS"]) < abs(density_only["U4:CDS"] - orf_data["truth"]["U4:CDS"])
