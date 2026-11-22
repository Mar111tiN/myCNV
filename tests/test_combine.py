"""Sample combine (step E) against main's combine_sample_CNV."""

import glob
import os

import numpy as np
import pytest

from conftest import RESTRICTED, RLEGACY, assert_frames_close, needs_restricted, read
from mycnv.combine import combine_sample, concat_tables
from mycnv.tools import natural_chroms


pytestmark = needs_restricted


@pytest.fixture
def pon_dir(tmp_path):
    os.symlink(f"{RLEGACY}/pon_mini_CNV", tmp_path / "CNV")
    return str(tmp_path)


@pytest.mark.parametrize("gc", [True, False])
def test_combine_matches_main_except_fullexonpos(pon_dir, gc):
    cov = concat_tables(sorted(glob.glob(f"{RESTRICTED}/sample_mini/*.cov.gz")))
    snp = concat_tables(sorted(glob.glob(f"{RESTRICTED}/sample_mini/*.snp")))
    new_cov, new_snp = combine_sample(cov, snp, pon_dir, gc_normalize=gc)
    tag = "gc" if gc else "nogc"
    old_cov, old_snp = read(f"{RLEGACY}/combine_{tag}.cov.gz"), read(f"{RLEGACY}/combine_{tag}.snp")
    assert_frames_close(new_cov, old_cov, exclude=("FullExonPos",))
    assert_frames_close(new_snp, old_snp, exclude=("FullExonPos",))
    # FullExonPos: offsets in natural chromosome order (main: lexicographic, chr21 before chr7)
    pon = read(f"{pon_dir}/CNV/pon.cov.filterGC.gz")
    chrom_max = pon.groupby("Chr")["ExonPos"].max()
    offset, add = 0, {}
    for chrom in natural_chroms(chrom_max.index):
        add[chrom] = offset
        offset += chrom_max[chrom]
    for df in (new_cov, new_snp):
        assert np.array_equal(df["FullExonPos"], df["ExonPos"] + df["Chr"].map(add))
    assert list(dict.fromkeys(new_cov["Chr"])) == ["chr7", "chr21", "chrX"]
