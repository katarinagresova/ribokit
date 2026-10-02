import filecmp

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from ribokit import annotation, quant
from ribokit.bam import _footprint_length

from synth import TRUE_OFFSETS, make_dataset

WINDOW = (28, 30)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("synth"))


@pytest.fixture(scope="module")
def run(data, tmp_path_factory):
    prefix = tmp_path_factory.mktemp("run") / "s"
    q, off, stats = quant.quantify(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(prefix))
    return dict(prefix=prefix, quant=q.set_index("Name"), offsets=off, stats=stats)


@pytest.mark.parametrize("which", ["gtf", "gtf_with_stop"])
def test_annotation_matches_truth(data, which):
    anno = annotation.load_annotation(str(data[which]), str(data["fasta"]))
    coding = {t: tx for t, tx in data["txs"].items() if tx["cds_start"] is not None}
    assert anno.tx == sorted(coding)
    for i, t in enumerate(anno.tx):
        tx = coding[t]
        assert (anno.tx_len[i], anno.cds_start[i], anno.cds_end[i]) == (len(tx["seq"]), tx["cds_start"], tx["cds_end"])
    assert anno.stats["gtf_cds_includes_stop"] == (which == "gtf_with_stop")
    assert anno.stats["no_cds"] == 1


def test_offsets_recovered(run):
    got = {(r.length, r.phase): r.offset for r in run["offsets"].itertuples()}
    assert got == TRUE_OFFSETS
    assert (run["offsets"]["z"] > 3).all()


def test_em_length_term():
    # design record A.2: A = 900 nt, B = 600 nt inside A, both 1 read/nt
    unique = np.array([300.0, 0.0])
    classes = sp.csr_matrix(np.array([[1.0, 1.0]]))
    alpha, *_ = quant.em(unique, classes, np.array([1200.0]), np.array([900.0, 600.0]), 1e-9, 100_000)
    np.testing.assert_allclose(alpha / alpha.sum(), [0.6, 0.4], atol=1e-6)


def test_counts_match_truth(run, data):
    q, truth = run["quant"], data["truth"]
    for t in ("A1", "D1"):
        assert abs(q.at[t, "NumReads"] - truth[t]) / truth[t] < 0.02
    b_share = q.at["B1", "NumReads"] / (q.at["B1", "NumReads"] + q.at["B2", "NumReads"])
    assert abs(b_share - truth["B1"] / (truth["B1"] + truth["B2"])) < 0.03
    assert q.at["C1", "NumReads"] == pytest.approx(q.at["C2", "NumReads"])
    assert abs(q.at["C1", "NumReads"] + q.at["C2", "NumReads"] - truth["C1"] - truth["C2"]) / (truth["C1"] + truth["C2"]) < 0.02
    assert q["NumReads"].sum() == pytest.approx(run["stats"]["reads_assigned"])
    assert q["ritpm"].sum() == pytest.approx(1e6)


def test_ties_listed(run):
    ties = pd.read_csv(f"{run['prefix']}.ties.tsv", sep="\t")
    assert sorted(ties["Name"]) == ["C1", "C2"]


def test_filters_counted(run, data):
    s = run["stats"]
    assert s["reads_forward"] == len(data["reads"])                   # one count per name, whatever its _x<n>
    assert s["alignments_dropped_cigar"] == 7                          # the 10M1D17M reads
    assert s["reads_cigar_ok"] == len(data["reads"]) - 7
    long_reads = sum(len(seq) == 35 or (cig.startswith("1S") and len(seq) == 36) for _, seq, cig, _ in data["reads"])
    assert s["reads_in_length_window"] == s["reads_cigar_ok"] - long_reads
    assert s["reads_on_cds_transcript"] == s["reads_in_length_window"] - 11  # the reads on "unannotated"


def test_footprint_length():
    assert _footprint_length([(4, 1), (0, 28)]) == 28      # 5' soft clip: not part of the footprint
    assert _footprint_length([(0, 28), (4, 2)]) == 30      # 3' soft clip: part of it
    assert _footprint_length([(0, 10), (2, 1), (0, 17)]) is None
    assert _footprint_length([(4, 3)]) is None


def test_deterministic_and_offsets_reusable(run, data, tmp_path):
    again = tmp_path / "again"
    quant.quantify(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(again))
    reused = tmp_path / "reused"
    quant.quantify(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(reused),
                   offsets_path=f"{run['prefix']}.offsets.tsv")
    for suffix in ("quant.tsv", "offsets.tsv", "ties.tsv", "stats.tsv"):
        assert filecmp.cmp(f"{run['prefix']}.{suffix}", f"{again}.{suffix}", shallow=False), suffix
    assert filecmp.cmp(f"{run['prefix']}.quant.tsv", f"{reused}.quant.tsv", shallow=False)
