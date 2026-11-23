"""Coverage sex calls (mycnv.sex): the detection, the PON chrX doubling and pon.json (M4), sex.tsv (M5)."""

import json
import logging
import os
import shutil

import numpy as np
import pandas as pd
import pytest

from conftest import LEGACY, SYNTH, assert_frames_close, read
from mycnv import sex
from mycnv.combine import (
    amazonize, combine_pon_coverage, combine_sample, concat_tables, legacy_male_columns, make_pon_coverage, normalize,
    normalize_gc, pon_coverage_and_sex,
)
from mycnv.config import default_config, load_config
from mycnv.pipeline import pon_combine, sex_step

CHROMS = ["chr7", "chr21", "chrX"]
PON = f"{SYNTH}/pon"
T = sex.thresholds()
PAR1, PAR2 = sex.PAR["chrX"]


def coverage_table(x_level=1.0, n_auto=600, n_x=300, n_par=0, par_level=1.0, seed=1, log2=False):
    """Bins on chr1..chr3 (level 1), chrX non-PAR (x_level) and PAR1 (par_level); coverage ~100 or log2ratio ~0."""
    rng = np.random.default_rng(seed)
    chrom = np.repeat(["chr1", "chr2", "chr3"], n_auto // 3).tolist() + ["chrX"] * (n_x + n_par)
    pos = list(100_000 + 10 * np.arange(n_auto)) + list(5_000_000 + 10 * np.arange(n_x)) \
        + list(PAR1[0] + 10 * np.arange(n_par))
    level = np.r_[np.ones(n_auto), np.full(n_x, x_level), np.full(n_par, par_level)]
    if log2:
        values = np.log2(level) + rng.normal(0, 0.15, len(level))
    else:
        values = 100 * level * rng.uniform(0.8, 1.2, len(level))
    return pd.DataFrame({"Chr": chrom, "Pos": pos, "v": values})


def one_call(df, log2=False, t=T):
    return sex.call_columns(df, ["v"], t, log2=log2)[0]


# ---------- the detection ------------------------------------------------------------------------------------------
@pytest.mark.parametrize("log2", [False, True])
def test_xx_and_xy(log2):
    xx, xy = one_call(coverage_table(1.0, log2=log2), log2), one_call(coverage_table(0.5, log2=log2), log2)
    assert xx["call"] == "XX" and abs(xx["chrx_ratio"] - 1) < 0.05
    assert xy["call"] == "XY" and abs(xy["chrx_ratio"] - 0.5) < 0.05 and abs(xy["log2_diff"] + 1) < 0.1
    assert xy["n_chrx_bins"] == 300 and xy["n_auto_bins"] == 600


def test_par_bins_are_excluded():
    # XY: two copies on the PAR; 400 PAR bins would move the chrX median to the XX level
    df = coverage_table(0.5, n_par=400, par_level=1.0)
    r = one_call(df)
    assert r["call"] == "XY" and r["n_chrx_bins"] == 300
    everything_x = df.assign(Pos=np.where(df["Chr"] == "chrX", 5_000_000, df["Pos"]))
    assert one_call(everything_x)["chrx_ratio"] > 0.85  # without the PAR exclusion it would look XX
    only_par = coverage_table(0.5, n_x=0, n_par=400)
    assert one_call(only_par)["call"] == "unknown" and one_call(only_par)["n_chrx_bins"] == 0


def test_par_boundaries_are_inclusive():
    pos = [PAR1[0] - 1, PAR1[0], PAR1[1], PAR1[1] + 1, PAR2[0] - 1, PAR2[0], PAR2[1], PAR2[1] + 1]
    df = pd.DataFrame({"Chr": ["chrX"] * len(pos) + ["chr1", "chr22", "chrY", "chr7_alt", "chrM"],
                       "Pos": pos + [PAR1[0]] * 5})
    x, auto = sex.region_masks(df)
    assert x.tolist() == [True, False, False, True, True, False, False, True] + [False] * 5
    assert auto.tolist() == [False] * 8 + [True, True, False, False, False]


def test_ambiguous_and_unknown():
    assert one_call(coverage_table(0.75))["call"] == "ambiguous"
    few = one_call(coverage_table(0.5, n_x=T["min_bins"] - 1))
    assert few["call"] == "unknown" and np.isfinite(few["chrx_ratio"])  # the ratio is still reported
    assert one_call(coverage_table(0.5, n_x=T["min_bins"]))["call"] == "XY"
    no_auto = coverage_table(0.5).query("Chr == 'chrX'")
    assert one_call(no_auto)["call"] == "unknown" and one_call(no_auto)["n_auto_bins"] == 0
    no_x = coverage_table(0.5, n_x=0)
    assert one_call(no_x)["call"] == "unknown" and np.isnan(one_call(no_x)["chrx_ratio"])


@pytest.mark.parametrize("log2", [False, True])
def test_nan_and_inf_are_ignored(log2):
    df = coverage_table(0.5, log2=log2)
    bad = df.index[df["Chr"] == "chrX"][:150]
    df.loc[bad[:100], "v"] = np.nan
    df.loc[bad[100:], "v"] = np.inf if not log2 else -np.inf
    r = one_call(df, log2)
    assert r["call"] == "XY" and r["n_chrx_bins"] == 150


def test_thresholds_from_config():
    cfg = load_config(overrides={"sex": {"xy_max": 0.8}})
    assert cfg["sex"] == {"xy_max": 0.8, "xx_min": 0.85, "min_bins": 100}
    t = sex.thresholds(cfg["sex"])
    assert one_call(coverage_table(0.75), t=t)["call"] == "XY"
    with pytest.raises(ValueError, match="xy_max <= xx_min"):
        sex.thresholds({"xy_max": 0.9, "xx_min": 0.85})
    with pytest.raises(ValueError, match="min_bins"):
        sex.thresholds({"min_bins": 0})


def test_config_defaults():
    assert default_config()["sex"] == {"xy_max": 0.65, "xx_min": 0.85, "min_bins": 100} == sex.DEFAULTS


def test_ascat_gender(caplog):
    row = {"sample": "N", "role": "normal", "chrx_ratio": 0.7, "n_chrx_bins": 500, "n_auto_bins": 9}
    tumor = dict(row, sample="T", role="tumor", call="XY")
    assert sex.ascat_gender([dict(row, call="XY"), dict(tumor, call="XX")]) == "XY"  # pairs: the normal's call
    assert sex.ascat_gender([tumor]) == "XY"  # no normal row (tumor-only): the tumor's call
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        assert sex.ascat_gender([dict(row, call="ambiguous"), tumor]) == "XX"
        assert sex.ascat_gender([dict(row, call="unknown"), tumor]) == "XX"
    assert sum("ASCAT gender XX" in r.message for r in caplog.records) == 2


# ---------- M4: PON columns ---------------------------------------------------------------------------------------
def test_pon_calls_double_the_same_columns_as_main():
    tables, calls = pon_coverage_and_sex(PON, chroms=CHROMS)
    by = {c["column"]: c for c in calls}
    assert {c for c, r in by.items() if r["call"] == "XY"} == {"Cov3", "Cov6"}
    assert {c for c, r in by.items() if r["call"] == "XX"} == {f"Cov{i}" for i in (1, 2, 4, 5, 7, 8)}
    assert all(abs(by[c]["chrx_ratio"] - 0.5) < 0.05 for c in ("Cov3", "Cov6"))
    assert {c for c, r in by.items() if r["doubled"]} == {"Cov3", "Cov6"}
    # main's rule on the same tables doubles the same columns: the coverage tables are unchanged
    cov = combine_pon_coverage(PON, CHROMS)
    assert legacy_male_columns(normalize_gc(cov)) == legacy_male_columns(normalize(cov)) == ["Cov3", "Cov6"]
    cov_cols = [f"Cov{i}" for i in range(1, 9)]
    for name, old in [("fullGC", amazonize(normalize_gc(cov))), ("full", amazonize(normalize(cov)))]:
        new = tables[name].set_index(["Chr", "Pos"])[cov_cols]
        old = old.astype({"Chr": str}).set_index(["Chr", "Pos"]).loc[new.index, cov_cols]
        assert_frames_close(new.reset_index(), old.reset_index())
    # the call is made before the doubling: the doubled tables would call every column XX
    after = sex.call_columns(tables["full"], cov_cols, T)
    assert {c["call"] for c in after} == {"XX"}


def _pon_copy(tmp_path, scale=None):
    """The synthetic PON's cov/ + snp/ in tmp_path/pon, optionally with chrX of some columns scaled."""
    pon = tmp_path / "pon"
    shutil.copytree(PON, pon)
    if scale:
        path = pon / "cov" / "chrX.cov.gz"
        df = read(path)
        for col, factor in scale.items():
            df[col] = df[col] * factor
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    return pon


def test_pon_json_records_the_calls(tmp_path):
    pon = _pon_copy(tmp_path)
    pon_combine(str(pon), str(tmp_path / "out"), chroms=CHROMS)
    meta = json.loads((tmp_path / "out" / "CNV" / "pon.json").read_text())
    assert [e["column"] for e in meta["sex"]] == [f"Cov{i}" for i in range(1, 9)]
    assert all(set(e) == {"column", "name", "chrx_ratio", "call", "doubled"} for e in meta["sex"])
    assert all(e["name"] is None for e in meta["sex"])  # no PON_list.txt
    assert [e["column"] for e in meta["sex"] if e["doubled"]] == ["Cov3", "Cov6"]
    assert [e["call"] for e in meta["sex"]] == ["XX", "XX", "XY", "XX", "XX", "XY", "XX", "XX"]
    assert all(round(e["chrx_ratio"], 4) == e["chrx_ratio"] for e in meta["sex"])
    th = meta["sex_thresholds"]
    assert {k: th[k] for k in ("xy_max", "xx_min", "min_bins")} == sex.DEFAULTS and th["doubling"] == "XY call"
    assert th["par_grch38"] == {"chrX": [[10_001, 2_781_479], [155_701_383, 156_030_895]]}


def test_pon_json_names_from_the_pon_list(tmp_path, caplog):
    pon = _pon_copy(tmp_path)
    (pon / "PON_list.txt").write_text("".join(f"bams/N{i}.bam\n" for i in range(1, 9)))
    pon_combine(str(pon), str(tmp_path / "a"), chroms=CHROMS)
    meta = json.loads((tmp_path / "a" / "CNV" / "pon.json").read_text())
    assert [e["name"] for e in meta["sex"]] == [f"bams/N{i}.bam" for i in range(1, 9)]
    (pon / "PON_list.txt").write_text("".join(f"N{i}.bam\n" for i in range(1, 8)))  # 7 entries, 8 columns
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        pon_combine(str(pon), str(tmp_path / "b"), chroms=CHROMS)
    meta = json.loads((tmp_path / "b" / "CNV" / "pon.json").read_text())
    assert all(e["name"] is None for e in meta["sex"]) and "PON_list.txt has 7 entries" in caplog.text


def test_ambiguous_and_unknown_columns_are_not_doubled(tmp_path, caplog):
    """Cov2 at 66 % chrX: ambiguous (not doubled, warning); main's rule (mean < 75) doubles it."""
    pon = _pon_copy(tmp_path, scale={"Cov2": 0.66})
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        tables, calls = pon_coverage_and_sex(str(pon), chroms=CHROMS)
    by = {c["column"]: c for c in calls}
    assert by["Cov2"]["call"] == "ambiguous" and not by["Cov2"]["doubled"]
    assert "PON Cov2: sex call ambiguous" in caplog.text
    cov = combine_pon_coverage(str(pon), CHROMS)
    x = cov["Chr"] == "chrX"
    norm = normalize(cov)
    full = tables["full"]
    assert np.allclose(full.loc[full["Chr"] == "chrX", "Cov2"], norm.loc[x, "Cov2"])  # not doubled
    assert "Cov2" in legacy_male_columns(norm)
    # --legacy-stats keeps main's rule (and says where it differs from the coverage call)
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        legacy_tables, legacy_calls = pon_coverage_and_sex(str(pon), chroms=CHROMS, legacy=True)
    assert {c["column"] for c in legacy_calls if c["doubled"]} == {"Cov2", "Cov3", "Cov6"}
    assert "main's chrX rule (--legacy-stats) doubles Cov2, Cov3, Cov6" in caplog.text
    lx = legacy_tables["full"].loc[legacy_tables["full"]["Chr"] == "chrX"]
    assert (lx["Pos"].to_numpy() == norm.loc[x, "Pos"].to_numpy()).all()
    kept = lx["Cov2"].notna().to_numpy()  # main's full tables carry the outlier NaNs (doubled Cov2 is high: many)
    assert kept.sum() > 100
    assert np.allclose(lx["Cov2"].to_numpy()[kept], 2 * norm.loc[x, "Cov2"].to_numpy()[kept])
    # a PON without chrX: every column unknown, nothing doubled, tables as before
    no_x = tmp_path / "nox"
    shutil.copytree(PON, no_x)
    os.remove(no_x / "cov" / "chrX.cov.gz")
    os.remove(no_x / "snp" / "chrX.snp.gz")
    _, calls = pon_coverage_and_sex(str(no_x), chroms=["chr7", "chr21"])
    assert {c["call"] for c in calls} == {"unknown"} and not any(c["doubled"] for c in calls)


def test_legacy_stats_tables_unchanged():
    """--legacy-stats still reproduces main's PON tables (the 75 rule), with the calls only recorded."""
    tables, calls = pon_coverage_and_sex(PON, chroms=CHROMS, legacy=True)
    for name in ["full", "fullGC", "filter", "filterGC"]:
        assert_frames_close(tables[name], read(f"{LEGACY}/pon_synth_CNV/pon.cov.{name}.gz"))
    assert {c["column"] for c in calls if c["doubled"]} == {"Cov3", "Cov6"}
    assert make_pon_coverage(PON, chroms=CHROMS, legacy=True).keys() == tables.keys()


# ---------- M5: samples -------------------------------------------------------------------------------------------
@pytest.fixture
def pon_dir(tmp_path):
    os.symlink(f"{LEGACY}/pon_synth_CNV", tmp_path / "CNV")
    return str(tmp_path)


def _raw(pon_dir, normal_x=1.0, tumor_x=1.0):
    files = sorted(os.listdir(f"{SYNTH}/sample"))
    cov = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".cov.gz")])
    snp = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".snp")])
    x = cov["Chr"] == "chrX"
    cov.loc[x, "Cov1"] *= normal_x
    cov.loc[x, "Cov2"] *= tumor_x
    return combine_sample(cov, snp, pon_dir)[0]


def test_sample_sex_table(pon_dir, tmp_path):
    out = sex_step(_raw(pon_dir), "T", "N", default_config(), str(tmp_path / "sex.tsv"))
    assert list(out.columns) == sex.SEX_COLUMNS
    assert out[["sample", "role", "column"]].values.tolist() == [["N", "normal", "log2ratio1"],
                                                                 ["T", "tumor", "log2ratio2"]]
    assert out["call"].tolist() == ["XX", "XX"] and out["ascat_gender"].tolist() == ["XX", "XX"]
    assert (out["n_chrx_bins"] > 400).all() and (out["n_auto_bins"] > 2000).all()
    assert (out["log2_diff"].abs() < 0.2).all()
    pd.testing.assert_frame_equal(read(tmp_path / "sex.tsv"), out)
    # XY normal (chrX halved): log2 difference ~ -1 against the effectively female PON; ASCAT gender XY
    xy = sex.pair_sex_table(_raw(pon_dir, normal_x=0.5, tumor_x=0.5), "T", "N")
    assert xy["call"].tolist() == ["XY", "XY"] and (xy["ascat_gender"] == "XY").all()
    assert ((xy["log2_diff"] + 1).abs() < 0.2).all() and ((xy["chrx_ratio"] - 0.5).abs() < 0.1).all()
    # a tumor that lost one X: its own call is XY, the pair's ASCAT gender stays the normal's XX
    loss = sex.pair_sex_table(_raw(pon_dir, tumor_x=0.5), "T", "N")
    assert loss["call"].tolist() == ["XX", "XY"] and (loss["ascat_gender"] == "XX").all()


def test_sample_without_chrx_is_unknown_and_xx(pon_dir, caplog):
    raw = _raw(pon_dir)
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        out = sex.pair_sex_table(raw.loc[raw["Chr"] != "chrX"], "T", "N")
    assert out["call"].tolist() == ["unknown", "unknown"] and (out["ascat_gender"] == "XX").all()
    assert (out["n_chrx_bins"] == 0).all() and out["chrx_ratio"].isna().all()
    assert "Sex call of the normal N is unknown" in caplog.text


def test_sex_table_needs_log2ratio_columns(pon_dir):
    raw = _raw(pon_dir).drop(columns="log2ratio2")
    with pytest.raises(ValueError, match="no log2ratio2 column"):
        sex.pair_sex_table(raw, "T", "N")
