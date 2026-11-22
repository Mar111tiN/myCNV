"""Public parity tier for steps C-H: the new code against main's outputs on the synthetic CNV PON + sample.

Inputs: tests/data/synth/{matrix,pon,sample} (tests/synth.py). Legacy outputs: tests/data/legacy/*synth*
(tests/legacy_outputs.py, unmodified main in the somvar-legacy env). These tests run without the restricted fixtures,
so CI and `docker build --target test` check values, not only shapes.
"""

import copy
import os
import subprocess
from io import StringIO

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, LEGACY, SYNTH, assert_frames_close, read
from mycnv.combine import combine_sample, compute_stats, concat_tables, make_pon_coverage, make_pon_snp, remove_outliers
from mycnv.convert import write_ascat
from mycnv.pipeline import rolling_step
from mycnv.tools import mawk, natural_chroms

CHROMS = ["chr7", "chr21", "chrX"]
PON = f"{SYNTH}/pon"
STATS = ["PONcov_mean", "PONcov_median", "PONcov_std"]
# columns interpolated from the coverage bins onto the SNP rows (mergeCov2SNP)
INTERPOLATED = ["log2ratio1", "log2ratio1_mean", "log2ratio2", "log2ratio2_mean", "PONcov_mean", "PONcov_median",
                "PONcov_std", "GCratio"]


@pytest.fixture
def legacy_config(study_config):
    """The study config with the exponent main actually used (it never read diff_exp, so always 2)."""
    cfg = copy.deepcopy(study_config)
    for section in cfg["rolling"].values():
        section["diff_exp"] = 2
    return cfg


# ---------- C: PON2CNV ------------------------------------------------------------------------------------------
def test_pon2cnv_legacy_mode_matches_main_up_to_the_first_base_fix(tmp_path):
    snp = tmp_path / "pon.snp"
    bed = f"{DATA}/haev7_covered_chr7.bed"
    filt = f"gunzip -c {SYNTH}/matrix/chr7.pon.gz | {mawk('filterBed')} {bed} -x -L -c chr7"
    cov = subprocess.run(["bash", "-o", "pipefail", "-c", f"{filt} | {mawk('PON2CNV')} -x -L -o {snp} -v 0.05 -d 50 -c 0"],
                         check=True, capture_output=True, text=True).stdout
    new, old = pd.read_csv(StringIO(cov), sep="\t"), read(f"{LEGACY}/pon2cnv_synth_chr7.cov")
    with open(snp) as a, open(f"{LEGACY}/pon2cnv_synth_chr7.snp") as b:
        assert a.read() == b.read()  # pooled SNPs: byte-identical
    assert (new[["Chr", "Pos", "ExonPos"]].to_numpy() == old[["Chr", "Pos", "ExonPos"]].to_numpy()).all()
    stream = pd.read_csv(StringIO(subprocess.run(["bash", "-o", "pipefail", "-c", filt], check=True,
                                                 capture_output=True, text=True).stdout), sep="\t")
    depth = stream["Depth"].str.split("=", expand=True)
    d = depth[0].str.split("|", expand=True).astype(int) + depth[1].str.split("|", expand=True).astype(int)
    keep = d.sum(axis=1) > 0
    stream, d = stream.loc[keep].reset_index(drop=True), d.loc[keep].reset_index(drop=True)
    first = (stream["Start"] // 10).drop_duplicates(keep="first").index
    expected = d.loc[first].to_numpy() / 10
    expected[0, :] = 0  # main kept the very first base of the stream
    cov_cols = [c for c in new.columns if c.startswith("Cov")]
    assert np.allclose(new[cov_cols].to_numpy() - old[cov_cols].to_numpy(), expected, atol=1e-9)


# ---------- D: PON combine --------------------------------------------------------------------------------------
def test_pon_combine_legacy_matches_main():
    tables = make_pon_coverage(PON, std_factor=2, chroms=CHROMS, legacy=True)
    for name in ["full", "fullGC", "filter", "filterGC"]:
        assert_frames_close(tables[name], read(f"{LEGACY}/pon_synth_CNV/pon.cov.{name}.gz"))
    assert_frames_close(make_pon_snp(PON, chroms=CHROMS), read(f"{LEGACY}/pon_synth_CNV/pon.snp.gz"))


def test_pon_combine_recomputes_filter_stats_and_lifts_male_x():
    new = make_pon_coverage(PON, std_factor=2, chroms=CHROMS)
    old = {n: read(f"{LEGACY}/pon_synth_CNV/pon.cov.{n}.gz") for n in ["full", "filter"]}
    assert_frames_close(new["full"][STATS], old["full"][STATS])
    cov_cols = [c for c in new["full"].columns if c.startswith("Cov")]
    assert new["full"][cov_cols].isna().sum().sum() == 0 < old["full"][cov_cols].isna().sum().sum()
    expected = compute_stats(remove_outliers(new["full"], std_factor=2))
    assert_frames_close(new["filter"][STATS], expected[STATS])
    assert not np.allclose(new["filter"]["PONcov_std"], old["filter"]["PONcov_std"], equal_nan=True)
    # the 2 male normals (Cov3, Cov6: half chrX coverage in the input) are lifted to the XX level
    from mycnv.combine import combine_pon_coverage, normalize_gc

    norm = normalize_gc(combine_pon_coverage(PON, CHROMS))
    before = norm.loc[norm["Chr"] == "chrX", cov_cols].mean()
    after = new["fullGC"].loc[new["fullGC"]["Chr"] == "chrX", cov_cols].mean()
    assert set(before[before < 75].index) == {"Cov3", "Cov6"}
    assert (after > 75).all() and np.allclose(after[["Cov3", "Cov6"]], 2 * before[["Cov3", "Cov6"]])


# ---------- E: sample combine -----------------------------------------------------------------------------------
@pytest.fixture
def pon_dir(tmp_path):
    os.symlink(f"{LEGACY}/pon_synth_CNV", tmp_path / "CNV")
    return str(tmp_path)


def _sample():
    files = sorted(os.listdir(f"{SYNTH}/sample"))
    cov = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".cov.gz")])
    snp = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".snp")])
    return cov, snp


@pytest.mark.parametrize("legacy_coords", [True, False])
@pytest.mark.parametrize("gc", [True, False])
def test_combine_matches_main(pon_dir, gc, legacy_coords):
    cov, snp = _sample()
    new_cov, new_snp = combine_sample(cov, snp, pon_dir, gc_normalize=gc, legacy_coords=legacy_coords)
    tag = "gc" if gc else "nogc"
    old_cov, old_snp = read(f"{LEGACY}/combine_synth_{tag}.cov.gz"), read(f"{LEGACY}/combine_synth_{tag}.snp")
    # main's natural chromosome order is lexicographic (chr21, chr7, chrX): compare in genomic order
    key = ["Chr", "Pos"]
    old_cov = old_cov.sort_values(key, key=lambda s: s.map({c: i for i, c in enumerate(CHROMS)}) if s.name == "Chr"
                                  else s).reset_index(drop=True)
    old_snp = old_snp.sort_values(key, key=lambda s: s.map({c: i for i, c in enumerate(CHROMS)}) if s.name == "Chr"
                                  else s).reset_index(drop=True)
    assert_frames_close(new_cov, old_cov, exclude=("FullExonPos",))
    # main merged the PON SNPs on the float map* columns too: SNPs whose map values differ from the PON's lost PONVAF
    lost = old_snp["PONVAF"].isna() & new_snp["PONVAF"].notna()
    assert 0 < lost.sum() < 0.1 * len(new_snp)
    assert_frames_close(new_snp.loc[~lost], old_snp.loc[~lost], exclude=("FullExonPos",))
    assert_frames_close(new_snp.loc[lost], old_snp.loc[lost], exclude=("FullExonPos", "PONVAF", "PONDepth"))
    # FullExonPos: offsets in natural chromosome order (main: lexicographic, chr21 before chr7)
    pon = read(f"{pon_dir}/CNV/pon.cov.filterGC.gz")
    chrom_max = pon.groupby("Chr")["ExonPos"].max()
    offset, add = 0, {}
    for chrom in natural_chroms(chrom_max.index):
        add[chrom] = offset
        offset += chrom_max[chrom]
    for df in (new_cov, new_snp):
        assert np.array_equal(df["FullExonPos"], df["ExonPos"] + df["Chr"].map(add))
    assert list(dict.fromkeys(new_cov["Chr"])) == CHROMS


def test_gc_normalisation_ignores_chromosomes_without_pon(pon_dir):
    """main normalised over chr1-22, X: extra contigs (chrY) must not change the autosomal GC factors."""
    cov, snp = _sample()
    y = cov.loc[cov["Chr"] == "chr21"].assign(Chr="chrY", Cov1=0.0, Cov2=5.0)
    with_y = pd.concat([cov.astype({"Chr": str}), y], ignore_index=True)
    a, _ = combine_sample(cov, snp, pon_dir)
    b, _ = combine_sample(with_y, snp, pon_dir)
    assert_frames_close(a, b)


# ---------- F/G: rolling ----------------------------------------------------------------------------------------
def _raw(chroms):
    cov, snp = read(f"{LEGACY}/combine_synth_gc.cov.gz"), read(f"{LEGACY}/combine_synth_gc.snp")
    return cov.loc[cov["Chr"].isin(chroms)], snp.loc[snp["Chr"].isin(chroms)]


def test_rolling_single_chromosome_matches_main(legacy_config):
    snp, cnv, _ = rolling_step(*_raw(["chr7"]), legacy_config)
    assert len(snp) > 30  # the BAF / ASCAT path gets real SNPs
    assert_frames_close(snp, read(f"{LEGACY}/rolling_synth_chr7.cnv.snp.gz"))
    old = read(f"{LEGACY}/rolling_synth_chr7.cnv.gz")
    key = ["Chr", "Pos", "VAF1"]
    assert_frames_close(cnv.sort_values(key, kind="stable").reset_index(drop=True),
                        old.sort_values(key, kind="stable").reset_index(drop=True))


def test_diff_exp_is_read(study_config, legacy_config):
    snp4, _, _ = rolling_step(*_raw(["chr7"]), study_config)  # diff_exp: 4 from the study config
    snp2, _, _ = rolling_step(*_raw(["chr7"]), legacy_config)
    assert_frames_close(snp4, snp2, exclude=("SNPdensityDiff",))
    assert np.allclose(snp4["SNPdensityDiff"], snp2["SNPdensityDiff"] ** 2, equal_nan=True)
    assert not np.allclose(snp4["SNPdensityDiff"], snp2["SNPdensityDiff"], equal_nan=True)


def test_rolling_three_chromosomes_documented_deltas(legacy_config):
    """Per-chromosome windows (CHANGELOG): main's SNPs are kept, at most a few more pass; the non-rolled columns are
    unchanged; SNPdensity/fallSNP move only by the cross-chromosome windows main had in its min-max normalisation."""
    new, _, _ = rolling_step(*_raw(CHROMS), legacy_config)
    old = read(f"{LEGACY}/rolling_synth_all.cnv.snp.gz")
    m = new.merge(old, on=["Chr", "Pos"], suffixes=("", "_old"), indicator=True, how="outer")
    assert (m["_merge"] != "right_only").all()  # every SNP main kept is kept
    assert (m["_merge"] == "left_only").sum() <= max(2, 0.05 * len(old))
    shared = m.loc[m["_merge"] == "both"]
    rolled = set(INTERPOLATED) | {"SNPdensity", "SNPdensityDiff", "fallSNP", "covLLH1", "covLLH2", "covLLH1_sum",
                                  "covLLH2_sum", "covLLH1_sumDiff", "covLLH2_sumDiff", "log2ratio1_meanDiff",
                                  "log2ratio2_meanDiff"}
    for col in [c for c in new.columns if c not in rolled and c not in ("Chr", "Pos")]:
        a, b = shared[col], shared[col + "_old"]
        if pd.api.types.is_numeric_dtype(a):
            assert np.allclose(a, b, rtol=1e-9, atol=1e-12, equal_nan=True), col
    for col in ["SNPdensity", "fallSNP"]:
        assert (shared[col] - shared[col + "_old"]).abs().max() < 0.01, col


# ---------- H: ASCAT input --------------------------------------------------------------------------------------
def test_write_ascat_matches_main(tmp_path):
    files = write_ascat(read(f"{LEGACY}/rolling_synth_chr7.cnv.snp.gz"), "syn_T-N", str(tmp_path))
    names = sorted(os.path.basename(f) for f in files)
    assert names == sorted(os.listdir(f"{LEGACY}/ascat_synth_chr7"))
    for path in files:
        legacy = os.path.join(LEGACY, "ascat_synth_chr7", os.path.basename(path))
        assert_frames_close(read(path), read(legacy), rtol=1e-12)
        with open(path) as a, open(legacy) as b:
            assert a.readline() == b.readline()  # header: "\tchrs\tpos\t<sample>"
    tumor = read(tmp_path / "syn_T-N_T_baf.tsv")
    assert tumor["chrs"].astype(str).eq("7").all() and not tumor["syn_T-N"].eq(0.5).all()
    assert read(tmp_path / "syn_T-N_germline_baf.tsv")["syn_T-N"].eq(0.5).all()
    assert read(tmp_path / "syn_T-N_germline_logr.tsv")["syn_T-N"].eq(0).all()
