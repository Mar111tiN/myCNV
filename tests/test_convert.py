"""ASCAT input (step H) against main's write_ASCAT."""

import os

import pytest

from conftest import RLEGACY, assert_frames_close, needs_restricted, read
from mycnv.convert import split_sample_name, write_ascat


@needs_restricted
def test_write_ascat_matches_main(tmp_path):
    files = write_ascat(read(f"{RLEGACY}/rolling_chr7.cnv.snp.gz"), "NHL1_A-B", str(tmp_path))
    names = sorted(os.path.basename(f) for f in files)
    assert names == sorted(os.listdir(f"{RLEGACY}/ascat_chr7"))
    for path in files:
        legacy = os.path.join(RLEGACY, "ascat_chr7", os.path.basename(path))
        # numeric comparison: main wrote floats that pandas 1.5 had parsed up to 1 ulp off
        assert_frames_close(read(path), read(legacy), rtol=1e-12)
        with open(path) as a, open(legacy) as b:
            assert a.readline() == b.readline()  # header: "\tchrs\tpos\t<sample>"


def test_sample_names():
    assert split_sample_name("NHL1_A-B") == ("A", "B")
    assert split_sample_name("IMMU_NSCLC_0001_T1-N1") == ("T1", "N1")  # main took the 2nd "_" field and failed
    with pytest.raises(ValueError):
        split_sample_name("sample")
