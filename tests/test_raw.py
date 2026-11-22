"""tn2cnv / pon2cnv (step B/C) against main's TN2CNV / PON2CNV."""

import numpy as np

from conftest import DATA, GC, GENMAP, LEGACY, RESTRICTED, RLEGACY, SYNTH, assert_frames_close, needs_restricted, read
from mycnv.raw import pon2cnv, tn2cnv


def test_tn2cnv_legacy_coords_match_main():
    cov, snp = tn2cnv(f"{SYNTH}/chr7.cpile", f"{DATA}/haev7_covered_chr7.bed", gc=GC, genmap=GENMAP, chrom="chr7",
                      legacy_coords=True)
    assert_frames_close(cov, read(f"{LEGACY}/tn2cnv_chr7.cov.gz"))
    assert_frames_close(snp, read(f"{LEGACY}/tn2cnv_chr7.snp"))


def test_tn2cnv_whole_sample_equals_chrom_run():
    a = tn2cnv(f"{SYNTH}/chr7.cpile", f"{DATA}/haev7_covered_chr7.bed", gc=GC, genmap=GENMAP)
    b = tn2cnv(f"{SYNTH}/chr7.cpile", f"{DATA}/haev7_covered_chr7.bed", gc=GC, genmap=GENMAP, chrom="chr7")
    assert_frames_close(a[0], b[0])
    assert_frames_close(a[1], b[1])


def test_tn2cnv_fixed_coords_keep_values():
    cov, snp = tn2cnv(f"{SYNTH}/chr7.cpile", f"{DATA}/haev7_covered_chr7.bed", gc=GC, genmap=GENMAP)
    old_cov, old_snp = read(f"{LEGACY}/tn2cnv_chr7.cov.gz"), read(f"{LEGACY}/tn2cnv_chr7.snp")
    assert len(cov) == len(old_cov)
    assert np.allclose(cov[["Cov1", "Cov2"]], old_cov[["Cov1", "Cov2"]])
    assert (np.diff(cov["ExonPos"]) > 0).all() and (np.diff(cov["Pos"]) > 0).all()
    assert_frames_close(snp, old_snp, exclude=("ExonPos",))


@needs_restricted
def test_pon2cnv_legacy_coords_match_main_up_to_first_base_fix():
    cov, snp = pon2cnv(f"{RESTRICTED}/chr7.pon.gz", f"{DATA}/haev7_padded_chr7.bed", gc=GC, genmap=GENMAP, chrom="chr7",
                       legacy_coords=True)
    old = read(f"{RLEGACY}/pon2cnv_py_chr7.cov.gz")
    cov_cols = [c for c in old.columns if c.startswith("Cov")]
    assert_frames_close(cov, old, exclude=cov_cols)
    assert (cov[cov_cols].to_numpy() >= old[cov_cols].to_numpy() - 1e-9).all()  # only additions (first base)
    assert_frames_close(snp, read(f"{RLEGACY}/pon2cnv_py_chr7.snp.gz"))
