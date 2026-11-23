"""End-to-end CLI: synthetic PONmatrix -> CNV PON -> per-pair chain (public fixtures only), option vocabulary,
failure modes (non-zero exit instead of empty outputs)."""

import os
import shutil
import stat
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import DATA, GC, GENMAP, SGC, SGENMAP, STUDY_CONFIG, SYNTH, assert_frames_close, read

MYCNV = [sys.executable, "-m", "mycnv.cli"]
BED7 = f"{DATA}/haev7_covered_chr7.bed"
TRACKS7 = ["--gc", GC, "--genmap", GENMAP]
MULTI_BED = f"{SYNTH}/multi.bed"
TRACKS = ["--gc", SGC, "--genmap", SGENMAP]
MULTI_MATRIX = [f"{SYNTH}/matrix_multi/{c}.pon.gz" for c in ["chr7", "chr21", "chrX"]]


def run(*args, check=True):
    res = subprocess.run([*MYCNV, *map(str, args)], capture_output=True, text=True)
    if check and res.returncode:
        raise AssertionError(f"mycnv {' '.join(map(str, args))} exited {res.returncode}:\n{res.stderr[-3000:]}")
    return res


def lines(path):
    return len(read(path)) if os.path.exists(path) else None


def test_version_and_config():
    assert "mycnv" in run("--version").stdout
    out = run("config").stdout
    assert "diff_exp: 4" in out and "GCnormalize: true" in out
    assert "pileup:" not in out and "ascat:" not in out  # not read by mycnv


def test_config_warns_about_unused_keys(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("CNV:\n  PONcoverage:\n    stdfactor: 3\n  ascat:\n    gender: XY\n")
    res = run("config", "--config", path)
    assert "PONcoverage.stdfactor" in res.stderr and "ascat" in res.stderr
    assert "stdfactor" not in res.stdout


def test_pon_and_sample_chain(tmp_path):
    pon = tmp_path / "PONmatrix"
    run("pon-coverage", "--ponmatrix", SYNTH, "--bed", BED7, *TRACKS7, "-o", pon)
    assert (pon / "cov" / "chr7.cov.gz").is_file() and (pon / "snp" / "chr7.snp.gz").is_file()
    run("pon-combine", "--ponmatrix", pon, "-o", pon, "--chrom", "chr7")
    for name in ["full", "fullGC", "filter", "filterGC"]:
        assert (pon / "CNV" / f"pon.cov.{name}.gz").is_file()
    out = tmp_path / "pair"
    res = run("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED7, *TRACKS7, "--ponmatrix", pon,
              "--sample", "S1_T-N", "-o", out, "--plots", check=False)
    assert res.returncode == 2 and "no SNP passed the SNP filter" in res.stderr  # no silent empty ASCAT input
    assert lines(out / "data" / "S1_T-N.cnv.gz") > 100 and not (out / "pre").exists()
    run("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED7, *TRACKS7, "--ponmatrix", pon,
        "--sample", "S1_T-N", "-o", out, "--plots", "--allow-empty")
    raw = read(out / "raw" / "S1_T-N.cov.gz")
    assert len(raw) > 150 and {"log2ratio1", "log2ratio2", "PONcov_mean", "FullExonPos"} <= set(raw.columns)
    assert raw["log2ratio1"].notna().mean() > 0.9  # sample and PON bins share coordinates
    cnv = read(out / "data" / "S1_T-N.cnv.gz")
    assert {"log2ratio1_mean", "log2ratio2_mean"} <= set(cnv.columns)
    for f in ["N_baf", "T_baf", "germline_baf", "N_logr", "T_logr", "germline_logr"]:
        assert (out / "pre" / f"S1_T-N_{f}.tsv").is_file()
    assert len(list((out / "plots").glob("*.jpg"))) >= 4


# ---------- the multi-chromosome synthetic kit: PON, split path, threads ------------------------------------------
@pytest.fixture(scope="module")
def multi_pon(tmp_path_factory):
    pon = tmp_path_factory.mktemp("multi") / "PON"
    run("pon-coverage", "--matrix", *MULTI_MATRIX, "--bed", MULTI_BED, *TRACKS, "-o", pon, "-t", 2)
    run("pon-combine", "--ponmatrix", pon, "-o", pon, "--chrom", "chr7", "--chrom", "chr21", "--chrom", "chrX")
    return pon


def test_pon_coverage_threads_give_the_same_tables(tmp_path, multi_pon):
    one = tmp_path / "one"
    # --matrix repeated and --chrom filtering: both forms accumulate
    run("pon-coverage", "--matrix", MULTI_MATRIX[0], "--matrix", *MULTI_MATRIX[1:], "--bed", MULTI_BED, *TRACKS,
        "-o", one, "--threads", 1)
    for chrom in ["chr7", "chr21", "chrX"]:
        for sub, ext in [("cov", "cov.gz"), ("snp", "snp.gz")]:
            pd.testing.assert_frame_equal(read(one / sub / f"{chrom}.{ext}"), read(multi_pon / sub / f"{chrom}.{ext}"))


def test_split_path_equals_sample(tmp_path, multi_pon):
    """tn2cnv per chromosome -> combine (all files) -> rolling -> ascat-input == `sample` (plan D16)."""
    whole = tmp_path / "whole"
    run("sample", "--pileup", f"{SYNTH}/multi.cpile", "--bed", MULTI_BED, *TRACKS, "--ponmatrix", multi_pon,
        "--sample", "M_T-N", "-o", whole, "--allow-empty")
    covs, snps = [], []
    for chrom in ["chr7", "chr21", "chrX"]:
        covs.append(tmp_path / f"{chrom}.cov.gz")
        snps.append(tmp_path / f"{chrom}.snp")
        run("tn2cnv", "--pileup", f"{SYNTH}/multi.cpile", "--bed", MULTI_BED, *TRACKS, "--chrom", chrom,
            "--out-cov", covs[-1], "--out-snp", snps[-1])
    with pytest.raises(AssertionError, match="no bins on 2 of 3 PON chromosomes"):
        run("combine", "--cov", covs[0], "--snp", snps[0], "--ponmatrix", multi_pon, "--out-cov", tmp_path / "x",
            "--out-snp", tmp_path / "y")
    run("combine", "--cov", covs[0], "--cov", *covs[1:], "--snp", *snps, "--ponmatrix", multi_pon,
        "--out-cov", tmp_path / "raw.cov.gz", "--out-snp", tmp_path / "raw.snp")
    run("rolling", "--cov", tmp_path / "raw.cov.gz", "--snp", tmp_path / "raw.snp", "--out-snp",
        tmp_path / "d.cnv.snp.gz", "--out-cnv", tmp_path / "d.cnv.gz")
    run("ascat-input", "--cnv", tmp_path / "d.cnv.snp.gz", "--sample", "M_T-N", "-o", tmp_path / "pre", "--allow-empty")
    pairs = [("raw.cov.gz", "raw/M_T-N.cov.gz"), ("raw.snp", "raw/M_T-N.snp"), ("d.cnv.gz", "data/M_T-N.cnv.gz"),
             ("d.cnv.snp.gz", "data/M_T-N.cnv.snp.gz")]
    for split, ref in pairs:
        pd.testing.assert_frame_equal(read(tmp_path / split), read(whole / ref))
    assert sorted(os.listdir(tmp_path / "pre")) == sorted(os.listdir(whole / "pre"))
    whole_cov = read(whole / "raw" / "M_T-N.cov.gz")
    assert list(dict.fromkeys(whole_cov["Chr"])) == ["chr7", "chr21", "chrX"]


def test_raw_mpileup_route(tmp_path):
    """--raw-pileup (cleanpileup -d added) == cleanpileup output given as a clean pileup; indels become I/i/D/d."""
    from mycnv.tools import mawk

    clean = tmp_path / "clean.cpile"
    with open(clean, "w") as fh:
        subprocess.run(["bash", "-o", "pipefail", "-c", f"{mawk('cleanpileup')} -d < {SYNTH}/raw.mpileup"],
                       stdout=fh, check=True)
    c = read(clean)
    assert list(c.columns) == ["Chr", "Start", "Ref", "Depth1", "Read1", "Depth2", "Read2"]
    assert c["Read1"].tolist() == [".,.,.,", "Ii..,.", ".d.*,.", "..Aa"]
    assert c["Read2"].tolist() == [".,..,", ".D.,.", "T,tT.", "Ii,."]
    for name, src, extra in [("raw", f"{SYNTH}/raw.mpileup", ["--raw-pileup"]), ("clean", clean, [])]:
        run("tn2cnv", "--pileup", src, "--bed", MULTI_BED, *TRACKS, "--out-cov", tmp_path / f"{name}.cov",
            "--out-snp", tmp_path / f"{name}.snp", "--min-depth", 4, *extra)
    pd.testing.assert_frame_equal(read(tmp_path / "raw.cov"), read(tmp_path / "clean.cov"))
    pd.testing.assert_frame_equal(read(tmp_path / "raw.snp"), read(tmp_path / "clean.snp"))
    assert len(read(tmp_path / "raw.snp")) == 2  # chr7:1002 (I/i) and chr21:501 (A/a), normal VAF 1/3 and 1/2


# ---------- option vocabulary ------------------------------------------------------------------------------------
def test_repeated_list_options_accumulate():
    from mycnv.cli import build_parser

    p = build_parser()
    a = p.parse_args(["pon-coverage", "--ponmatrix", "m", "--bed", "b", "--gc", "g", "-o", "o", "--chrom", "chr1",
                      "--chrom", "chr2"])
    assert a.chrom == ["chr1", "chr2"]
    a = p.parse_args(["combine", "--cov", "a", "--cov", "b", "c", "--snp", "s", "--snp", "t", "--ponmatrix", "p",
                      "--out-cov", "x", "--out-snp", "y"])
    assert a.cov == ["a", "b", "c"] and a.snp == ["s", "t"]
    a = p.parse_args(["pon-combine", "--ponmatrix", "p", "-o", "o", "--chrom", "chr7", "--chrom", "chr7"])
    from mycnv.cli import _check

    _check(p, a)
    assert a.chrom == ["chr7"]


@pytest.mark.parametrize("argv", [
    ["sample", "--tumor", "T"],                   # ebscore's --tumor is a BAM; no prefix matching of --tumor-name
    ["pon-combine", "--ponmatrix", "p", "-o", "o", "--chroms", "chr7"],
    ["sample", "--pileup", "p", "--ponmatrix", "d", "--sample", "s", "-o", "o", "--bed", "b", "--genmap", "m"],
    ["pon-combine", "--ponmatrix", "p"],          # -o/--outdir is required: inputs are never written to
])
def test_rejected_arguments(argv):
    res = run(*argv, check=False)
    assert res.returncode == 2 and "error" in res.stderr


def test_threads_and_quiet_on_every_command():
    from mycnv.cli import build_parser

    p = build_parser()
    sub_type = type(p._subparsers._group_actions[0])

    def leaves(parser, prefix=""):
        subs = [a for a in parser._actions if isinstance(a, sub_type)]
        if not subs:
            yield prefix.strip(), parser
            return
        for name, child in subs[0].choices.items():
            yield from leaves(child, f"{prefix} {name}")

    names = []
    for name, parser in leaves(p):
        names.append(name)
        opts = {o for a in parser._actions for o in a.option_strings}
        assert {"-t", "--threads", "-v", "--verbose", "--quiet", "--config"} <= opts, name
    assert {"tracks gc", "tracks genmap", "sample"} <= set(names)
    sample = dict(leaves(p))["sample"]
    assert "--chrom" not in {o for a in sample._actions for o in a.option_strings}


def test_plots_dir_needs_sample(tmp_path):
    res = run("rolling", "--cov", "x", "--snp", "y", "--out-snp", tmp_path / "a", "--out-cnv", tmp_path / "b",
              "--plots-dir", tmp_path / "p", check=False)
    assert res.returncode == 2 and "--plots-dir needs --sample" in res.stderr


def test_one_name_alone_is_rejected(tmp_path):
    res = run("ascat-input", "--cnv", "x", "--sample", "s", "-o", tmp_path, "--tumor-name", "T", check=False)
    assert res.returncode == 2 and "--normal-name" in res.stderr


# ---------- failure modes: non-zero exit, nothing silently empty ---------------------------------------------------
@pytest.fixture(scope="module")
def pon7(tmp_path_factory):
    pon = tmp_path_factory.mktemp("pon7")
    run("pon-coverage", "--ponmatrix", SYNTH, "--bed", BED7, *TRACKS7, "-o", pon)
    run("pon-combine", "--ponmatrix", pon, "-o", pon, "--chrom", "chr7")
    return pon


def _sample(pon, out, *extra, pileup=f"{SYNTH}/chr7.cpile", tracks=TRACKS7):
    return run("sample", "--pileup", pileup, "--bed", BED7, *tracks, "--ponmatrix", pon, "--sample", "S1_T-N",
               "-o", out, "--allow-empty", *extra, check=False)


@pytest.mark.parametrize("tracks, message", [
    (["--gc", f"{DATA}/tracks/gc-typo", "--genmap", GENMAP], "does not exist"),
    (["--gc", GC, "--genmap", f"{DATA}/tracks/genmap-typo"], "does not exist"),
    (["--gc", SGC, "--genmap", GENMAP], "no GC data for chr7|no value at any"),
])
def test_wrong_tracks_fail(tmp_path, pon7, tracks, message):
    import re

    res = _sample(pon7, tmp_path / "out", tracks=tracks)
    assert res.returncode == 2 and re.search(message, res.stderr)
    assert not (tmp_path / "out" / "raw").exists()


def test_pon_coverage_with_wrong_gc_fails(tmp_path):
    res = run("pon-coverage", "--ponmatrix", SYNTH, "--bed", BED7, "--gc", f"{DATA}/tracks/GC", "-o", tmp_path / "p",
              check=False)
    assert res.returncode == 2 and "does not exist" in res.stderr


def test_header_only_pileup_fails(tmp_path, pon7):
    empty = tmp_path / "empty.cpile"
    empty.write_text("Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2\n")
    res = _sample(pon7, tmp_path / "out", pileup=empty)
    assert res.returncode == 2 and "no coverage bins" in res.stderr


def test_read_only_ponmatrix_and_paths_with_spaces(tmp_path, pon7):
    """Staged inputs are read-only and keep their (possibly odd) names; outputs go to -o only."""
    staged = tmp_path / "sp ace" / "PON matrix"
    shutil.copytree(pon7, staged)
    pile = tmp_path / "sp ace" / "pair 1.cpile"
    shutil.copy(f"{SYNTH}/chr7.cpile", pile)
    for root, dirs, files in os.walk(staged):
        for name in dirs + files:
            os.chmod(os.path.join(root, name), stat.S_IREAD | (stat.S_IEXEC if name in dirs else 0))
    os.chmod(staged, stat.S_IREAD | stat.S_IEXEC)
    try:
        run("pon-combine", "--ponmatrix", staged, "-o", tmp_path / "cnv pon", "--chrom", "chr7")
        res = _sample(tmp_path / "cnv pon", tmp_path / "out dir", pileup=pile)
        assert res.returncode == 0, res.stderr[-2000:]
        assert lines(tmp_path / "out dir" / "raw" / "S1_T-N.cov.gz") > 150
    finally:
        for root, dirs, files in os.walk(staged):
            os.chmod(root, stat.S_IRWXU)
            for name in files:
                os.chmod(os.path.join(root, name), stat.S_IRUSR | stat.S_IWUSR)


# ---------- the packaged defaults are the study parameters and reach the steps ------------------------------------
def test_defaults_equal_the_study_config_on_every_key_read():
    from mycnv.config import default_config

    with open(STUDY_CONFIG) as fh:
        study = yaml.safe_load(fh)["CNV"]
    d = default_config()
    keys = [("hetSNP", "normalVAF"), ("hetSNP", "minDepth"), ("coverage", "minCov"), ("coverage", "GCnormalize"),
            ("PONcoverage", "minCov"), ("PONcoverage", "minVAF"), ("PONcoverage", "minDepth"),
            ("PONcoverage", "stdFactor"), ("filter", "cov"), ("filter", "snp")]
    for sec in ["cov", "snp"]:
        keys += [("rolling", sec, k) for k in ["normalize", "debug", "diff_exp", "ddof"]]
    keys += [("rolling", "cov", "LLH"), ("rolling", "cov", "data"), ("rolling", "snp", "data", "SNPdensity")]
    for path in keys:
        a, b = d, study
        for k in path:
            a, b = a[k], b[k]
        assert a == b, path


def test_options_reach_the_steps(tmp_path, pon7):
    """--normal-vaf / --min-depth change the het SNPs; --no-gc-normalize uses the non-GC PON table."""
    base = ["tn2cnv", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED7, *TRACKS7]
    run(*base, "--out-cov", tmp_path / "a.cov", "--out-snp", tmp_path / "a.snp")
    run(*base, "--out-cov", tmp_path / "b.cov", "--out-snp", tmp_path / "b.snp", "--normal-vaf", 0.4, 0.6)
    run(*base, "--out-cov", tmp_path / "c.cov", "--out-snp", tmp_path / "c.snp", "--min-depth", 100)
    a, b, c = (read(tmp_path / f"{x}.snp") for x in "abc")
    assert len(b) < len(a) and b["VAF1"].between(0.4, 0.6).all() and len(a) > 0
    assert len(c) < len(a) and (c[["Depth1", "Depth2"]] >= 100).all().all()
    for name, extra in [("gc", []), ("nogc", ["--no-gc-normalize"])]:
        _sample(pon7, tmp_path / name, *extra)
    gc, nogc = read(tmp_path / "gc" / "raw" / "S1_T-N.cov.gz"), read(tmp_path / "nogc" / "raw" / "S1_T-N.cov.gz")
    for df, table in [(gc, "filterGC"), (nogc, "filter")]:
        pon = read(pon7 / "CNV" / f"pon.cov.{table}.gz")
        m = df[["Pos", "PONcov_mean"]].merge(pon[["Pos", "PONcov_mean"]], on="Pos")
        assert np.allclose(m["PONcov_mean_x"], m["PONcov_mean_y"])
    assert not np.allclose(gc["PONcov_mean"], nogc["PONcov_mean"])
    assert not np.allclose(gc["Cov1"], nogc["Cov1"])
