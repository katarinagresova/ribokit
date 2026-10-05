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
    # window [c - 1, c + 3) inside the background, no start within 8 nt
    assert got == {c for c in range(65, 148) if abs(c - 100) > 8}


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


def test_orf_types():
    cs, ce = np.full(8, 100), np.full(8, 400)
    start = np.array([10, 40, 97, 98, 100, 160, 161, 420])
    end = np.array([60, 400, 400, 130, 400, 400, 200, 500])
    assert list(starts.orf_types(start, end, cs, ce)) == [
        "uORF", "extension", "extension", "uoORF", "CDS", "truncation", "internal", "dORF"]
    assert starts.orf_types(np.array([5]), np.array([50]), np.array([-1]), np.array([-1]))[0] == "other"


def test_bh_counts_untested_as_p1():
    p = np.array([0.001, 0.02, 0.5])
    assert np.allclose(starts.bh(p, 3), sps.false_discovery_control(p))
    assert np.allclose(starts.bh(p, 10), sps.false_discovery_control(np.r_[p, np.ones(7)])[:3])


def test_near_starts_best_first():
    pos = np.array([100, 106, 112, 150, 200])
    p = np.array([1e-3, 1e-8, 1e-5, 1e-6, 0.5])
    called, by = starts.near_starts(pos, p, p < 0.01)
    # 106 is best: it stops 100 (6 nt before) and 112 (codon 2 after); 150 is 44 nt after, out of reach
    assert called.tolist() == [False, True, False, True, False]
    assert by.tolist() == [1, -1, 1, -1, -1]
    # codons 1-10 of an annotated start (here 100) are not called, even when it is not called itself
    called, by = starts.near_starts(pos, p, p < 0.01, annotated=[0])
    assert called.tolist() == [True, False, False, True, False] and by.tolist() == [-1, 0, 0, -1, -1]


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
    scan, scan_st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(d / "scan"), scan=True)
    return dict(dir=d, data=data, sheet=sheet, out=out, stats=st, scan=scan, scan_stats=scan_st)


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
    cli.main(["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
              "--scan", "--out-prefix", str(tmp_path / "scan")])
    for f in ("starts", "start_orfs", "starts_stats"):
        assert filecmp.cmp(d / f"scan.{f}.tsv", tmp_path / f"scan.{f}.tsv", shallow=False), f


def test_scan_calls(started):
    # check 1 for the scan. Over seeds 3-12: every annotated start called; 61 of 70 planted starts called at their
    # codon, the misses at 4-5x enrichment with 1-5 elongation reads (q 0.06-0.3); the pair's CTG and the weaker
    # start 5 codons after a planted one are never called, and listed with the start that stops them; the pause
    # and the untranslated uORF never called. Calls at no planted start: 0-5 per seed in leaders (q >= 0.012 but
    # one), 2-8 in CDSs (about BH's share: the CDS null windows are at 0.6-0.8x alpha)
    data = started["data"]
    out = started["scan"]
    key = {f"{n}_{s + 1}": i for i, (n, s, e, _, _) in data["orf_coords"].items()}
    out = out.assign(truth=out["ORF_id"].map(key)).set_index("ORF_id")
    by_truth = out.dropna(subset=["truth"]).set_index("truth")
    cds = by_truth[by_truth.index.str.endswith(":CDS")]
    assert cds["called"].all() and (cds["type"] == "CDS").all()
    assert sum(by_truth.at[i, "called"] for i in PLANTED) >= 6
    want = {"S1_uATG": "uORF", "S1_uCTG": "uORF", "S2_uoORF": "uoORF", "S3_ext": "extension", "S4_pairATG": "uORF",
            "S5_first": "uORF", "N1_ORF": "other"}
    for i, t in want.items():
        assert by_truth.at[i, "type"] == t, i
        assert by_truth.at[i, "end"] == data["orf_coords"][i][2], i
    for weaker, stronger in (("S4_pairCTG", "S4_pairATG"), ("S5_after", "S5_first")):
        if weaker in by_truth.index:
            assert not by_truth.at[weaker, "called"]
            if by_truth.at[stronger, "called"]:
                assert out.at[by_truth.at[weaker, "stopped_by"], "truth"] == stronger
    for i in ("S6_pause", "S7_quiet"):
        assert i not in by_truth.index or not by_truth.at[i, "called"], i
    false_leader = out[out["called"] & out["truth"].isna() & (out["region"] == "leader")]
    assert len(false_leader) <= 3
    so = pd.read_csv(started["dir"] / "scan.start_orfs.tsv", sep="\t")
    assert len(so) == started["scan_stats"]["called"] - len(cds) - started["scan_stats"]["called_without_stop"]
    assert (so["type"] != "CDS").all() and set(so.columns[:4]) == {"ORF_id", "Name", "start", "end"}


@pytest.mark.parametrize("table", ["start_orfs", "catalogue"])
def test_scan_orfs_feed_orfs(started, tmp_path, table):
    # start_orfs.tsv and catalogue.tsv are ORF tables for `ribokit orfs`: every row is kept, and the catalogue's CDS
    # rows are the annotated CDSs (none added)
    d, data = started["dir"], started["data"]
    _, _, _, st = orfs.run(str(data["libs"]["elongation", 1]), str(data["gtf"]), str(data["fasta"]), WINDOW,
                           str(tmp_path / "o"), orfs_path=str(d / f"scan.{table}.tsv"),
                           offsets_path=str(d / "pool.offsets.tsv"))
    assert st["orf_table_rows"] == len(pd.read_csv(d / f"scan.{table}.tsv", sep="\t"))
    assert sum(v for k, v in st.items() if k.startswith("orfs_dropped_")) == 0
    if table == "catalogue":
        assert st["orfs_cds_added"] == 0


def test_catalogue(started):
    out = started["scan"]
    cat = pd.read_csv(started["dir"] / "scan.catalogue.tsv", sep="\t").set_index("ORF_id")
    calls = out[out["called"] & out["type"].isin(["uORF", "uoORF"])]
    assert set(cat.index[cat["type"] != "CDS"]) == set(calls["ORF_id"])
    assert (cat["type"] == "CDS").sum() == 37
    data = started["data"]
    for i in ("S1_uATG", "S2_uoORF", "N1_ORF", "S3_ext"):
        n, s, e, _, _ = data["orf_coords"][i]
        assert (f"{n}_{s + 1}" in cat.index) == (i in ("S1_uATG", "S2_uoORF")), i   # no extension, no ORF without a CDS
