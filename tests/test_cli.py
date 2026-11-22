"""End-to-end CLI: synthetic PONmatrix -> CNV PON -> per-pair chain (public fixtures only)."""

import subprocess
import sys

import pandas as pd

from conftest import DATA, GC, GENMAP, SGC, SGENMAP, SYNTH

MYCNV = [sys.executable, "-m", "mycnv.cli"]


def run(*args):
    return subprocess.run([*MYCNV, *map(str, args)], check=True, capture_output=True, text=True)


def test_version_and_config():
    assert "mycnv" in run("--version").stdout
    out = run("config").stdout
    assert "diff_exp: 4" in out and "GCnormalize: true" in out


def test_pon_and_sample_chain(tmp_path):
    pon = tmp_path / "PONmatrix"
    bed = f"{DATA}/haev7_covered_chr7.bed"  # one BED flavour for PON and samples (plan L9)
    run("pon-coverage", "--matrix-dir", f"{SYNTH}/matrix", "--bed", bed, "--gc", GC, "--genmap", GENMAP,
        "--out-dir", pon)
    assert (pon / "cov" / "chr7.cov.gz").is_file() and (pon / "snp" / "chr7.snp.gz").is_file()
    run("pon-combine", "--pon-dir", pon, "--chroms", "chr7")
    for name in ["full", "fullGC", "filter", "filterGC"]:
        assert (pon / "CNV" / f"pon.cov.{name}.gz").is_file()
    out = tmp_path / "pair"
    run("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", bed, "--gc", GC, "--genmap", GENMAP, "--pon-dir", pon,
        "--sample", "S1_T-N", "--out-dir", out)
    raw = pd.read_csv(out / "raw" / "S1_T-N.cov.gz", sep="\t")
    assert len(raw) > 150 and {"log2ratio1", "log2ratio2", "PONcov_mean", "FullExonPos"} <= set(raw.columns)
    assert raw["log2ratio1"].notna().mean() > 0.9  # sample and PON bins share coordinates
    cnv = pd.read_csv(out / "data" / "S1_T-N.cnv.gz", sep="\t")
    assert {"log2ratio1_mean", "log2ratio2_mean"} <= set(cnv.columns)
    for f in ["N_baf", "T_baf", "germline_baf", "N_logr", "T_logr", "germline_logr"]:
        assert (out / "pre" / f"S1_T-N_{f}.tsv").is_file()


def test_split_runs_equal_whole_sample(tmp_path):
    """tn2cnv per chromosome + combine == the whole-sample run (optional split, plan D16)."""
    bed = f"{SYNTH}/multi.bed"
    tracks = ["--gc", SGC, "--genmap", SGENMAP]
    run("tn2cnv", "--pileup", f"{SYNTH}/multi.cpile", "--bed", bed, *tracks, "--out-cov", tmp_path / "all.cov",
        "--out-snp", tmp_path / "all.snp")
    parts = []
    for chrom in ["chr7", "chr21", "chrX"]:
        run("tn2cnv", "--pileup", f"{SYNTH}/multi.cpile", "--bed", bed, *tracks, "--chrom", chrom,
            "--out-cov", tmp_path / f"{chrom}.cov", "--out-snp", tmp_path / f"{chrom}.snp")
        parts.append(pd.read_csv(tmp_path / f"{chrom}.cov", sep="\t"))
    whole = pd.read_csv(tmp_path / "all.cov", sep="\t")
    pd.testing.assert_frame_equal(whole, pd.concat(parts, ignore_index=True))
