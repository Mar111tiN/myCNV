"""PON combine (step D) against main's make_PON_coverage / make_PON_snp."""

import numpy as np

import pytest

from conftest import RESTRICTED, RLEGACY, assert_frames_close, needs_restricted, read
from mycnv.combine import (
    combine_pon_coverage, compute_stats, legacy_male_columns, make_pon_coverage, make_pon_snp, normalize, normalize_gc,
    pon_coverage_and_sex, remove_outliers,
)

PON = f"{RESTRICTED}/pon_mini"
pytestmark = needs_restricted
CHROMS = ["chr7", "chr21", "chrX"]


def test_pon_combine_legacy_matches_main():
    tables = make_pon_coverage(PON, std_factor=2, chroms=CHROMS, legacy=True)
    for name in ["full", "fullGC", "filter", "filterGC"]:
        assert_frames_close(tables[name], read(f"{RLEGACY}/pon_mini_CNV/pon.cov.{name}.gz"))
    assert_frames_close(make_pon_snp(PON, chroms=CHROMS), read(f"{RLEGACY}/pon_mini_CNV/pon.snp.gz"))


def test_pon_combine_recomputes_filter_stats():
    new = make_pon_coverage(PON, std_factor=2, chroms=CHROMS)
    old = {n: read(f"{RLEGACY}/pon_mini_CNV/pon.cov.{n}.gz") for n in ["full", "filter"]}
    stats = ["PONcov_mean", "PONcov_median", "PONcov_std"]
    # full tables: same unfiltered stats as main, but no outlier NaNs in the coverage
    assert_frames_close(new["full"][stats], old["full"][stats])
    cov_cols = [c for c in new["full"].columns if c.startswith("Cov")]
    assert new["full"][cov_cols].isna().sum().sum() < old["full"][cov_cols].isna().sum().sum()
    # filter tables: stats recomputed on the outlier-free coverage (main kept the unfiltered ones)
    expected = compute_stats(remove_outliers(new["full"], std_factor=2))
    assert_frames_close(new["filter"][stats], expected[stats])
    assert not np.allclose(new["filter"]["PONcov_std"], old["filter"]["PONcov_std"], equal_nan=True)


def test_sex_calls_double_the_same_columns_as_main():
    """pon_mini's chrX slice (chrX:284,050-1,465,540) lies wholly in PAR1, where XY normals have two copies: every
    column is unknown (0 non-PAR bins) and none is doubled; main's 75 rule doubles none either. (So this check is
    vacuous; the real-data comparison was run on the whole HAEv7 PON, see CHANGELOG "Sex calls measured".)"""
    tables, calls = pon_coverage_and_sex(PON, chroms=CHROMS)
    cov = combine_pon_coverage(PON, CHROMS)
    assert legacy_male_columns(normalize_gc(cov)) == legacy_male_columns(normalize(cov)) == []
    assert len(calls) == 28 and {c["call"] for c in calls} == {"unknown"} and {c["n_chrx_bins"] for c in calls} == {0}
    assert not any(c["doubled"] for c in calls)
