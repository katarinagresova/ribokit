import filecmp
import itertools

import numpy as np
import pandas as pd
import pytest

from ribokit import orfs, score

from synth import make_orf_dataset

WINDOW = (28, 30)


def test_lead_probability_exact():
    rng = np.random.default_rng(0)
    for n in range(1, 7):
        for q in [np.full(3, 1 / 3), rng.dirichlet(np.ones(3)), np.array([0.9, 0.1, 0.0])]:
            brute = sum(np.prod(q[list(fs)]) for fs in itertools.product(range(3), repeat=n)
                        if fs.count(0) > max(fs.count(1), fs.count(2)))
            assert score.lead_probability(n, q) == pytest.approx(brute, abs=1e-12), (n, q)


def test_upper_tail_exact():
    p = np.array([0.1, 0.5, 0.33, 0.9, 0.02])
    for k in range(7):
        brute = sum(np.prod(np.where(x, p, 1 - p)) for x in itertools.product([0, 1], repeat=p.size) if sum(x) >= k)
        assert score.upper_tail(p, k) == pytest.approx(brute, abs=1e-12), k


@pytest.fixture(scope="module")
def scored(tmp_path_factory):
    # background runs under the ORFs too, so the untranslated candidate has reads in random frames
    d = tmp_path_factory.mktemp("scoresynth")
    data = make_orf_dataset(d, background_under_orfs=True)
    orfs.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(d / "s"), orfs_path=str(data["orfs"]))
    scores, decoys, stats = score.run([str(d / "s")], str(d / "out"), frame_lengths=WINDOW)
    return dict(dir=d, scores=scores.set_index("ORF_id"), decoys=decoys, stats=stats)


def test_synthetic_calls(scored):
    # check 1: translated ORFs are called; the untranslated candidate, the CUG with a start peak only
    # and U4 (untranslated, but the CDS start peak spills into it) are not. Over seeds 2-11, the
    # untranslated three had q >= 0.088, the translated ones q <= 0.057
    q = scored["scores"]["q"]
    for i in ("U1_uORF", "U2_uoORF", "U3_uoORF", "N1_ORF"):
        assert q[i] < 0.05, i
    for i in ("U1_cand", "U1_CUG", "U4_uoORF"):
        assert q[i] > 0.05, i
    assert scored["scores"].at["U1_CUG", "codons"] == 10       # its start codon, the peak, has no vote
    # decoys: translated ORFs' are not called; U2 and U4 (frame 1) lose +2, U3 (frame 2) +1: in the CDS's frame
    d = scored["decoys"]
    assert (d.loc[d["ORF_id"].isin(["U1_uORF", "U2_uoORF", "U3_uoORF", "N1_ORF"]), "p"] > 0.05).all()
    assert set(zip(d["ORF_id"], d["shift"])).isdisjoint({("U2_uoORF", 2), ("U3_uoORF", 1), ("U4_uoORF", 2)})
    assert scored["stats"]["decoys_left_out_cds_frame"] == 3


def test_pooled_and_deterministic(scored):
    d = scored["dir"]
    score.run([str(d / "s")], str(d / "again"), frame_lengths=WINDOW)
    for suffix in ("scores.tsv", "decoys.tsv", "score_stats.tsv"):
        assert filecmp.cmp(d / f"out.{suffix}", d / f"again.{suffix}", shallow=False), suffix
    # a library pooled with itself: twice the reads, the same leads
    twice, _, _ = score.run([str(d / "s")] * 2, str(d / "twice"), frame_lengths=WINDOW)
    twice = twice.set_index("ORF_id")
    assert (twice["reads"] == 2 * scored["scores"]["reads"]).all()
    assert (twice["leads"] == scored["scores"]["leads"]).all()


def test_null_follows_host_cds():
    # an untranslated in-frame extension of a CDS: its 100 codons in the CDS lead, as the CDS's
    # reads do. The null mixes in the CDS's frame profile, so it is not called; with a flat 1/3
    # null (no CDS reads) it would be
    o = pd.DataFrame({"ORF_id": ["T:CDS", "T_ext"], "Name": "T", "type": ["CDS", "other"], "start": [90, 30],
                      "end": 390, "Length": [300, 360], "NumReads": [1000.0, 50.0]})
    frames = pd.DataFrame({"length": [28], "offset": [12], "reads": [1000.0], "frame0": [0.9], "frame1": [0.05],
                           "frame2": [0.05]})
    upstream = [(j, *np.eye(3, dtype=int)[j % 3]) for j in range(1, 20)]     # background, one read per codon
    in_cds = [(j, 9, 1, 0) for j in range(20, 120)]
    codons = pd.DataFrame(upstream + in_cds, columns=["codon", "frame0", "frame1", "frame2"]).assign(
        ORF_id="T_ext", length=28)
    with_cds, _, _ = score.score([(codons, o, frames)], [28])
    assert with_cds.at[0, "leads"] == 106 and with_cds.at[0, "p"] > 0.3
    flat, _, _ = score.score([(codons, o.assign(NumReads=[np.nan, 50.0]), frames)], [28])
    assert flat.at[0, "p"] < 1e-10
