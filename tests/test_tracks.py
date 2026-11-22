"""GC / mappability tracks fail loudly instead of giving empty outputs; per-chromosome tables concatenate safely."""

import logging
import os
import shutil

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, GC, GENMAP, SGC, SGENMAP, SYNTH, read
from mycnv.combine import concat_tables
from mycnv.raw import pon2cnv, tn2cnv
from mycnv.tracks import add_gc, add_genmap, load_gc

BED = f"{DATA}/haev7_covered_chr7.bed"
COV = pd.DataFrame({"Chr": "chr7", "Pos": np.arange(10000, 10200, 10), "ExonPos": np.arange(20) * 10 + 1})


def test_gc_track_is_required_and_must_exist(tmp_path):
    with pytest.raises(ValueError, match="no GC track"):
        add_gc(COV, None)
    with pytest.raises(FileNotFoundError, match="does not exist"):
        add_gc(COV, str(tmp_path / "gc-typo"))
    with pytest.raises(ValueError, match="no mappability track"):
        tn2cnv(f"{SYNTH}/chr7.cpile", BED, gc=GC, genmap=None)


def test_missing_chromosome_file_fails_for_primary_chromosomes_only(tmp_path, caplog):
    empty = tmp_path / "gc"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="no GC data for chr7"):
        add_gc(COV, str(empty))
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        out = add_gc(COV.assign(Chr="chrY"), str(empty))
    assert out["GCratio"].isna().all() and "no GC data for chrY" in caplog.text


def test_track_that_matches_no_bin_fails(tmp_path):
    gc_dir = tmp_path / "gc"
    gc_dir.mkdir()
    gc = read(f"{GC}/chr7.gc100-10.gz")
    gc["Start"] += 5  # off the 10 bp grid: no bin start matches
    gc.to_csv(gc_dir / "chr7.gc100-10.gz", sep="\t", index=False, compression="gzip")
    with pytest.raises(ValueError, match="no value at any of the 20 bins"):
        add_gc(COV, str(gc_dir))


def test_single_gc_file_without_chr_column_needs_a_matching_name(tmp_path):
    path = f"{GC}/chr7.gc100-10.gz"
    assert len(load_gc(path, "chr7")) > 0
    with pytest.raises(ValueError, match="does not name chr1"):
        load_gc(path, "chr1")
    copy = tmp_path / "gc.tsv.gz"
    shutil.copy(path, copy)
    with pytest.raises(ValueError, match="no Chr column"):
        add_gc(COV, str(copy))


def test_single_genmap_file_with_chr_column():
    gm = read(f"{SGENMAP}/hg38_genmap.chr21.txt.gz")
    one = os.path.join(SGENMAP, "hg38_genmap.chr21.txt.gz")
    df = pd.DataFrame({"Chr": "chr21", "Pos": gm["Pos"].iloc[:30].to_numpy()})
    (out,) = add_genmap(df, track=one)
    assert out["map30_0"].notna().all()
    with pytest.raises(FileNotFoundError, match="no mappability data for chr7"):
        add_genmap(df.assign(Chr="chr7"), track=one)


def test_pon_coverage_does_not_need_genmap():
    cov, snp = pon2cnv(f"{SYNTH}/matrix/chr7.pon.gz", BED, gc=GC, genmap=None)
    assert len(cov) and cov["map30_0"].isna().all() and cov["GCratio"].notna().all()


def test_empty_pileup_fails(tmp_path):
    path = tmp_path / "empty.cpile"
    path.write_text("Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2\n")
    with pytest.raises(ValueError, match="no coverage bins"):
        tn2cnv(str(path), BED, gc=GC, genmap=GENMAP)
    with pytest.raises(ValueError, match="no region on chr22"):
        tn2cnv(f"{SYNTH}/chr7.cpile", BED, gc=GC, genmap=GENMAP, chrom="chr22")


def test_multi_chromosome_tracks_and_chromosome_list():
    cov, snp = tn2cnv(f"{SYNTH}/multi.cpile", f"{SYNTH}/multi.bed", gc=SGC, genmap=SGENMAP, chrom=["chrX", "21", "chrX"])
    assert list(dict.fromkeys(cov["Chr"])) == ["chr21", "chrX"]
    assert cov["GCratio"].notna().all() and cov["map30_0"].notna().all()


def test_concat_keeps_numeric_dtypes_with_a_header_only_file(tmp_path):
    cov, _ = tn2cnv(f"{SYNTH}/multi.cpile", f"{SYNTH}/multi.bed", gc=SGC, genmap=SGENMAP, chrom="chr21")
    cov.to_csv(tmp_path / "a.cov", sep="\t", index=False)
    cov.iloc[:0].to_csv(tmp_path / "b.cov", sep="\t", index=False)
    df = concat_tables([tmp_path / "b.cov", tmp_path / "a.cov"])
    assert len(df) == len(cov)
    for col in ["Pos", "ExonPos", "GCratio", "Cov1", "Cov2"]:
        assert pd.api.types.is_numeric_dtype(df[col]), col
