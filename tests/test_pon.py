"""PON combine (step D) against main's make_PON_coverage / make_PON_snp."""

import numpy as np

import pytest

from conftest import RESTRICTED, RLEGACY, assert_frames_close, needs_restricted, read
from mycnv.combine import compute_stats, make_pon_coverage, make_pon_snp, remove_outliers

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
