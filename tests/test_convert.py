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


@pytest.mark.parametrize("name", ["P01-T_vs_P01-N", "T1-1_vs_N1-1", "T-01_vs_N-01", "P1_-N", "P1_T-"])
def test_pair_ids_that_are_not_legacy_names_are_rejected(name):
    with pytest.raises(ValueError, match="--tumor-name"):
        split_sample_name(name)


def _cnv(n=3, logr=0.1):
    import pandas as pd

    return pd.DataFrame({"Chr": ["chr7"] * n, "Pos": range(100, 100 + n), "VAF1": 0.5, "VAF2": 0.4,
                         "log2ratio1_mean": logr, "log2ratio2_mean": 0.2})


def test_one_name_alone_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="both"):
        write_ascat(_cnv(), "P1_T-N", str(tmp_path), tumor="TUMOR")


def test_explicit_names_and_roles(tmp_path):
    files = write_ascat(_cnv(), "P01-T_vs_P01-N", str(tmp_path), tumor="P01-T", normal="P01-N")
    names = sorted(os.path.basename(f) for f in files)
    assert "P01-T_vs_P01-N_P01-T_baf.tsv" in names and "P01-T_vs_P01-N_P01-N_logr.tsv" in names
    tumor = read(tmp_path / "P01-T_vs_P01-N_P01-T_baf.tsv")
    assert (tumor["P01-T_vs_P01-N"] == 0.4).all() and (tumor["chrs"].astype(str) == "7").all()


def test_empty_or_logr_less_ascat_input_fails_unless_allowed(tmp_path):
    with pytest.raises(ValueError, match="no SNP passed"):
        write_ascat(_cnv(0), "P1_T-N", str(tmp_path / "a"))
    assert not (tmp_path / "a").exists()
    with pytest.raises(ValueError, match="has no value"):
        write_ascat(_cnv(3, logr=float("nan")), "P1_T-N", str(tmp_path / "b"))
    files = write_ascat(_cnv(0), "P1_T-N", str(tmp_path / "c"), allow_empty=True)
    assert len(files) == 6
