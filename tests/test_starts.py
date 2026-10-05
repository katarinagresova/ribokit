import filecmp

import numpy as np
import pandas as pd
import pysam
import pytest
from scipy import stats as sps

from ribokit import cli, orfs, starts
from ribokit.quant import quantify

from synth import make_start_dataset

WINDOW = (28, 30)
PLANTED = ["S1_uATG", "S1_uCTG", "S2_uoORF", "S3_ext", "S4_pairATG", "S5_first", "N1_ORF"]


def test_nb_sf_matches_scipy():
    rng = np.random.default_rng(0)
    h = rng.integers(0, 30, 500)
    mean = rng.random(500) * 20
    size = rng.random(500) * 5 + 0.1
    size[:100] = np.inf
    fin = np.isfinite(size)
    ref = np.where(fin, sps.nbinom.sf(h - 1, np.where(fin, size, 1), np.where(fin, size, 1) / (np.where(fin, size, 1) + mean)),
                   sps.poisson.sf(h - 1, mean))
    assert np.allclose(starts.nb_sf(h, mean, size), ref, atol=1e-12)


def test_null_codons_keep_clear_of_starts():
    bg = np.zeros(200, dtype=bool)
    bg[64:150] = True
    is_start = np.zeros(200, dtype=bool)
    is_start[100] = True
    got = set(starts.null_codons(bg, is_start).tolist())
    # window [c - 3, c + 6) inside the background, no start within 8 nt
    assert got == {c for c in range(67, 145) if abs(c - 100) > 8}


def test_fit_factors_recover_planted():
    rng = np.random.default_rng(1)
    n_tx, per = 400, 60
    tx = np.repeat(np.arange(n_tx), per)
    b = np.tile(np.repeat(np.arange(3), per // 3), n_tx)
    lam = rng.gamma(2, 1, tx.size)
    f_true, s_true = np.array([4.0, 1.0, 0.5]), rng.uniform(0.5, 2, n_tx)
    e = rng.poisson(lam)
    h = rng.poisson(s_true[tx] * f_true[b] * lam)
    pos = np.arange(tx.size)
    f, s = starts.fit_factors(h, e, tx, b, pos, pos, n_tx, 3, prior_reads=0)
    assert f / f[1] == pytest.approx(f_true, rel=0.1)


def test_start_lengths_longest_run():
    peaks = pd.DataFrame({"a": [9, 0, 0, -3, 0, 0, 0], "b": [9, 0, 9, 0, 0, 0, -1]}, index=range(24, 31))
    assert starts.choose_start_lengths(peaks) == [28, 29]


@pytest.fixture(scope="module")
def started(tmp_path_factory):
    # as in the validation: one offset table per pool, from `quant` on the merged harringtonine BAMs
    d = tmp_path_factory.mktemp("startsynth")
    data = make_start_dataset(d)
    pysam.cat("-o", str(d / "harr.bam"), *[str(data["libs"]["harringtonine", r]) for r in (1, 2)])
    quantify(str(d / "harr.bam"), str(data["gtf"]), str(data["fasta"]), WINDOW, str(d / "pool"))
    rows = []
    for rep in (1, 2):
        for kind in ("harringtonine", "elongation"):
            orfs.run(str(data["libs"][kind, rep]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(d / f"{kind}{rep}"),
                     orfs_path=str(data["orfs"]), offsets_path=str(d / "pool.offsets.tsv"))
        rows.append(("A", str(d / f"harringtonine{rep}"), str(d / f"elongation{rep}")))
    sheet = d / "sheet.tsv"
    pd.DataFrame(rows, columns=["pool", "harringtonine", "elongation"]).to_csv(sheet, sep="\t", index=False)
    out, st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(d / "s"), str(data["orfs"]))
    return dict(dir=d, data=data, sheet=sheet, out=out, stats=st)


def test_start_lengths_from_kernel(started):
    assert started["stats"]["start_lengths"] == "28,29,30"
    k = pd.read_csv(started["dir"] / "s.kernel.tsv", sep="\t")
    peak = k.loc[k.groupby("length")["kernel"].idxmax(), "nt"]
    assert (peak == 0).all()


def test_synthetic_calls(started):
    # check 1 for given starts. Over seeds 3-12: every annotated start called; 2 of 70 planted starts missed
    # (q 0.07-0.09, a uORF ATG and the pair's ATG); the pause and the untranslated uORF never called; at most
    # one of the ten untranslated filler candidates called (q 0.009-0.025: BH among true starts). The pair's
    # CTG gets the ATG's nt -1 reads and is called in 6 of 10 seeds: the scan's near-start rule, not this
    # test, drops it. The pair's ATG and the first of the two close starts always have the smaller p.
    out = started["out"].set_index("ORF_id")
    q, p = out["q"], out["p"]
    cds = out[out.index.str.endswith(":CDS")]
    assert (cds["q"] < 0.05).all()
    assert sum(q[i] < 0.05 for i in PLANTED) >= 6
    assert q["S6_pause"] > 0.05 and q["S7_quiet"] > 0.05
    assert (q[out.index.str.endswith("_cand")] < 0.05).sum() <= 1
    assert p["S4_pairATG"] < p["S4_pairCTG"] and p["S5_first"] < p["S5_after"]
    assert out.at["S6_pause", "enrichment"] < 1.5        # a pause: elongation peak, harringtonine not raised
    assert out.at["S2_uoORF", "region"] == "leader" and pd.isna(out.at["N1_ORF", "frame"])


def test_starts_cli_deterministic(started, tmp_path):
    d = started["dir"]
    data = started["data"]
    cli.main(["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
              "--orfs", str(data["orfs"]), "--out-prefix", str(tmp_path / "again")])
    for f in ("starts", "kernel", "factors", "starts_stats"):
        assert filecmp.cmp(d / f"s.{f}.tsv", tmp_path / f"again.{f}.tsv", shallow=False), f
