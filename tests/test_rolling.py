"""Rolling step (filter, rolling coverage, rolling SNP, remerge) on NHL1 raw slices against main."""

import copy

import numpy as np
import pandas as pd
import pytest

from conftest import RESTRICTED, RLEGACY, assert_frames_close, needs_restricted, read
from mycnv.pipeline import rolling_step
from mycnv.rolling import one_col_rolling

# columns interpolated from the coverage bins onto the SNP rows (mergeCov2SNP)
INTERPOLATED = ["log2ratio1", "log2ratio1_mean", "log2ratio2", "log2ratio2_mean", "PONcov_mean", "PONcov_median",
                "PONcov_std", "GCratio"]


def _run(tag, cfg):
    return rolling_step(read(f"{RESTRICTED}/nhl1/raw_{tag}.cov.gz"), read(f"{RESTRICTED}/nhl1/raw_{tag}.snp"), cfg)


@pytest.fixture
def legacy_config(study_config):
    """The study config with the exponent main actually used (it never read diff_exp, so always 2)."""
    cfg = copy.deepcopy(study_config)
    for section in cfg["rolling"].values():
        section["diff_exp"] = 2
    return cfg


def _tie_adjacent(snp: pd.DataFrame, cov: pd.DataFrame) -> np.ndarray:
    """SNP rows whose nearest coverage bin (left or right in FullExonPos) shares its FullExonPos with another bin.

    The legacy ExonPos is not monotonic, so bins of different regions can share a FullExonPos; main's unstable
    sort then picked either bin as interpolation neighbour. The new code keeps genomic (Pos) order on ties.
    """
    pos = np.sort(cov["FullExonPos"].to_numpy())
    dup = set(cov.loc[cov["FullExonPos"].duplicated(keep=False), "FullExonPos"])
    i = np.searchsorted(pos, snp["FullExonPos"].to_numpy())
    left = pos[np.clip(i - 1, 0, len(pos) - 1)]
    right = pos[np.clip(i, 0, len(pos) - 1)]
    return np.array([l in dup or r in dup for l, r in zip(left, right)])


@needs_restricted
def test_rolling_single_chromosome_matches_main(legacy_config):
    snp, cnv, _ = _run("chr7", legacy_config)
    old_snp = read(f"{RLEGACY}/rolling_chr7.cnv.snp.gz")
    assert_frames_close(snp, old_snp, exclude=INTERPOLATED)
    ties = _tie_adjacent(snp, read(f"{RESTRICTED}/nhl1/raw_chr7.cov.gz"))
    for col in INTERPOLATED:
        same = np.isclose(snp[col], old_snp[col], rtol=1e-9, atol=1e-12, equal_nan=True)
        assert (same | ties).all(), f"{col} differs on a row without a FullExonPos tie"
        assert same[~ties].all()
    assert ties.sum() <= 3  # documented: 1 of 9 SNP rows on this slice
    # full track: same rows; coverage rows identical
    old_cnv = read(f"{RLEGACY}/rolling_chr7.cnv.gz")
    key = ["Chr", "Pos", "VAF1"]
    a = cnv.sort_values(key, kind="stable").reset_index(drop=True)
    b = old_cnv.sort_values(key, kind="stable").reset_index(drop=True)
    cov_rows = a["VAF1"].isna().to_numpy()
    assert_frames_close(a.loc[cov_rows], b.loc[cov_rows])
    assert_frames_close(a.loc[~cov_rows], b.loc[~cov_rows], exclude=INTERPOLATED)


@needs_restricted
def test_diff_exp_is_read(study_config, legacy_config):
    snp4, _, _ = _run("chr7", study_config)  # diff_exp: 4 from the study config
    snp2, _, _ = _run("chr7", legacy_config)
    # only the Diff column changes: Diff**4 == (Diff**2)**2
    assert_frames_close(snp4, snp2, exclude=("SNPdensityDiff",))
    assert np.allclose(snp4["SNPdensityDiff"], snp2["SNPdensityDiff"] ** 2, equal_nan=True)


def test_rolling_does_not_cross_chromosomes():
    rng = np.random.default_rng(0)
    n = 1000
    df = pd.DataFrame({
        "Chr": ["chr7"] * n + ["chr21"] * n,
        "Pos": np.r_[np.arange(n), np.arange(n)] * 10,
        "FullExonPos": np.arange(2 * n) * 10,
        "log2ratio1": np.r_[rng.normal(0, 0.05, n), rng.normal(1, 0.05, n)],
    })
    out = one_col_rolling(df, "log2ratio1", "mean", 250, {"normalize": False, "diff_exp": 2})
    a = out.loc[out["Chr"] == "chr7", "log2ratio1_mean"]
    b = out.loc[out["Chr"] == "chr21", "log2ratio1_mean"]
    # main's windows mixed both levels over ~250 rows at the boundary; per chromosome they stay apart
    assert a.max() < 0.2 and b.min() > 0.8


def test_rolling_short_chromosome_uses_its_length():
    df = pd.DataFrame({"Chr": ["chr21"] * 40, "Pos": np.arange(40) * 10, "FullExonPos": np.arange(40) * 10,
                       "log2ratio1": np.linspace(-0.1, 0.1, 40)})
    out = one_col_rolling(df, "log2ratio1", "mean", 250, {"normalize": False, "diff_exp": 2})
    assert out["log2ratio1_mean"].notna().all()


@needs_restricted
def test_rolling_chromosome_order_does_not_matter(legacy_config):
    cov, snp = read(f"{RESTRICTED}/nhl1/raw_chr7_21.cov.gz"), read(f"{RESTRICTED}/nhl1/raw_chr7_21.snp")
    a, _, _ = rolling_step(cov, snp, legacy_config)
    # shift chr21's FullExonPos far up (main's lexicographic offsets put chr21 before chr7)
    for df in (cov, snp):
        df.loc[df["Chr"] == "chr21", "FullExonPos"] += 10**8
    b, _, _ = rolling_step(cov, snp, legacy_config)
    assert_frames_close(a, b, exclude=("FullExonPos",))
