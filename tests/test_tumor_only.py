"""Tumor-only PON-CNV (design item M1, plan D27 / N-Q3, D-Q3 mode B): one-sample pileup, het-SNP candidates at given
common-SNP sites, one sample column downstream, mock germline, one-row sex.tsv. Public synthetic inputs only: the
one-sample pileup is the tumor columns of synth/multi.cpile, the sites file is made from it (tests/synth.py)."""

import gzip
import logging
import os
import re
import subprocess
import sys
from io import StringIO

import numpy as np
import pandas as pd
import pytest

import synth
from conftest import LEGACY, SGC, SGENMAP, SYNTH, assert_frames_close, expected_het_snps, read
from mycnv import sex
from mycnv.combine import combine_sample, concat_tables
from mycnv.config import default_config, load_config
from mycnv.convert import write_ascat
from mycnv.pipeline import check_sample_columns, cnv_plots, raw_plots, rolling_step, sex_step
from mycnv.raw import pileup_samples, tn2cnv, tumor2cnv
from mycnv.sites import check_sites, sites_at
from mycnv.tools import mawk

MYCNV = [sys.executable, "-m", "mycnv.cli"]
MULTI_BED = f"{SYNTH}/multi.bed"
MULTI_PILE = f"{SYNTH}/multi.cpile"
TRACKS = ["--gc", SGC, "--genmap", SGENMAP]
MULTI_MATRIX = [f"{SYNTH}/matrix_multi/{c}.pon.gz" for c in ["chr7", "chr21", "chrX"]]
HET_EVERY = 7  # synth.main: make_pileup(MULTI, seed=2, het_every=7)
BAND = (0.1, 0.9)  # hetSNP.tumorOnlyVAF default
SNP_COLS = ["Chr", "Pos", "ExonPos", "Depth1", "VAF1"]


def run(*args, check=True):
    res = subprocess.run([*MYCNV, *map(str, args)], capture_output=True, text=True)
    if check and res.returncode:
        raise AssertionError(f"mycnv {' '.join(map(str, args))} exited {res.returncode}:\n{res.stderr[-3000:]}")
    return res


def write_gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def stream(pileup):
    """The filterBed stream (with ExonPos) of a pileup over multi.bed: the rows pile2CNV sees."""
    out = subprocess.run(["bash", "-o", "pipefail", "-c", f"{mawk('filterBed')} {MULTI_BED} -x < {pileup}"],
                         check=True, capture_output=True, text=True).stdout
    return pd.read_csv(StringIO(out), sep="\t")


def keys(df):
    return set(zip(df["Chr"].astype(str), df["Pos"].astype(int)))


@pytest.fixture(scope="module")
def tumor(tmp_path_factory):
    """(folder, one-sample pileup, its text, site records, VCF.gz of the records) from synth/multi.cpile."""
    root = tmp_path_factory.mktemp("tumor_only")
    with open(MULTI_PILE) as fh:
        text = synth.make_tumor_pileup(fh.read())
    pile = root / "M_T.cpile"
    pile.write_text(text)
    records = synth.make_snp_sites(text, het_every=HET_EVERY)
    vcf = write_gz(root / "sites.vcf.gz", synth.sites_vcf(records))
    return root, pile, text, records, vcf


# ---------- het-SNP candidates at the sites (raw.tumor2cnv, mycnv.sites) -----------------------------------------
def test_tumor_only_het_snps_are_the_sites_within_band_and_depth(tmp_path, tumor):
    _, pile, text, records, _ = tumor
    st = stream(pile)
    cand = expected_het_snps(st, normal_vaf=BAND, min_depth=15)  # tumor = sample 1: the pile2CNV rule on the tumor
    c21 = cand.loc[cand["Chr"] == "chr21", "Pos"].astype(int).tolist()
    multi, indel, dup, split, split_indel = c21[:5]
    special = {("chr21", multi): ("G", "A,T"), ("chr21", indel): ("GA", "G")}  # not biallelic SNVs: dropped
    recs = [(c, p, *special.get((c, p), (r, a))) for c, p, r, a in records]
    recs.append(next(r for r in recs if r[:2] == ("chr21", dup)))  # an identical record listed twice: counted once
    # a multi-allelic site split into biallelic records (bcftools norm -m-, gnomAD v4) is dropped like a joined one:
    # a second SNV allele, or an indel at the same position
    ref = {r[:2]: r[2] for r in recs}
    recs.append(("chr21", split, ref["chr21", split], next(b for b in "CGT" if b != ref["chr21", split])))
    recs.append(("chr21", split_indel, ref["chr21", split_indel], ref["chr21", split_indel] + "TT"))
    dropped = set(special) | {("chr21", split), ("chr21", split_indel)}
    vcf = write_gz(tmp_path / "sites.vcf.gz", synth.sites_vcf(recs))
    cov, snp = tumor2cnv(str(pile), MULTI_BED, str(vcf), gc=SGC, genmap=SGENMAP)
    assert all(p % HET_EVERY == 0 for p in cand["Pos"])  # in-band tumor candidates are the het positions
    sites = {(c, p) for c, p, *_ in recs} - dropped
    exp = cand.loc[[k in sites for k in zip(cand["Chr"], cand["Pos"].astype(int))]]
    got = snp[SNP_COLS].astype({"Pos": int, "ExonPos": int, "Depth1": int, "VAF1": float})
    assert_frames_close(got, exp.reset_index(drop=True), rtol=1e-5)
    # what each rule removed
    assert len(snp) >= 20 and keys(snp) <= sites
    assert snp["VAF1"].between(*BAND).all() and (snp["Depth1"] >= 15).all()
    assert ((snp["VAF1"] < 0.2) | (snp["VAF1"] > 0.8)).any()  # the tumor band, not the pair's normal band
    assert len(cand) > len(snp)  # tumor candidates that are no sites
    in_pileup = sites & keys(st.rename(columns={"Start": "Pos"}))  # sites inside the BED
    hom = {k for k in in_pileup if k[1] % 101 == 0 and k[1] % HET_EVERY}  # tumor VAF 0.95
    assert len(hom) == 5 and not hom & keys(snp)
    ref_like = {k for k in in_pileup if k[1] % 101 and k[1] % HET_EVERY}  # tumor VAF < 0.1
    assert len(ref_like) > 10 and not ref_like & keys(snp)
    assert not dropped & keys(snp)
    assert int(((snp["Chr"] == "chr21") & (snp["Pos"] == dup)).sum()) == 1
    # the depth rule
    _, deep = tumor2cnv(str(pile), MULTI_BED, str(vcf), gc=SGC, genmap=SGENMAP, min_depth=60)
    exp60 = expected_het_snps(st, normal_vaf=BAND, min_depth=60)
    assert keys(deep) == keys(exp60) & sites and 0 < len(deep) < len(snp) and (deep["Depth1"] >= 60).all()
    # a narrower band
    _, narrow = tumor2cnv(str(pile), MULTI_BED, str(vcf), gc=SGC, genmap=SGENMAP, tumor_vaf=(0.4, 0.6))
    assert 0 < len(narrow) < len(snp) and narrow["VAF1"].between(0.4, 0.6).all()
    assert list(snp.columns) == ["Chr", "Pos", "ExonPos", "map30_0", "map50_0", "map75_1", "map100_2", "Depth1", "VAF1"]


def test_tumor_only_coverage_bins_are_the_pairs(tumor):
    """Coverage bins and ExonPos come from the same filterBed | pile2CNV pass: Cov1 == the pair's tumor Cov2."""
    _, pile, _, _, vcf = tumor
    cov_t, _ = tumor2cnv(str(pile), MULTI_BED, str(vcf), gc=SGC, genmap=SGENMAP)
    cov_p, _ = tn2cnv(MULTI_PILE, MULTI_BED, gc=SGC, genmap=SGENMAP)
    pd.testing.assert_frame_equal(cov_t, cov_p.drop(columns="Cov1").rename(columns={"Cov2": "Cov1"}))
    cov_c, _ = tumor2cnv(str(pile), MULTI_BED, str(vcf), gc=SGC, genmap=SGENMAP, chrom="chr21")
    pd.testing.assert_frame_equal(cov_c, cov_t.loc[cov_t["Chr"] == "chr21"].reset_index(drop=True))


def test_sites_formats_give_the_same_snps(tmp_path, tumor):
    """VCF plain/gzipped (also without the ##fileformat line) and a Chr Pos TSV with or without a header (`Chr Pos`,
    or a # header such as the `#CHROM POS` of `grep -v '^##' sites.vcf | cut -f1,2`) select the same sites."""
    _, pile, _, records, vcf = tumor
    plain = tmp_path / "sites.vcf"
    plain.write_text(synth.sites_vcf(records))
    bare = tmp_path / "bare.vcf"  # only the #CHROM POS ID REF ALT ... column header
    bare.write_text("".join(line for line in synth.sites_vcf(records).splitlines(True) if not line.startswith("##")))
    body = "".join(f"{c}\t{p}\n" for c, p, *_ in records)
    files = [vcf, plain, bare]
    for name, head in [("h.tsv", "Chr\tPos\n"), ("n.tsv", ""), ("c.tsv", "#Chr\tPos\n"), ("v.tsv", "#CHROM\tPOS\n")]:
        (tmp_path / name).write_text(head + body)
        files.append(tmp_path / name)
    files.append(write_gz(tmp_path / "h.tsv.gz", "Chr\tPos\n" + body))
    ref = None
    for path in files:
        _, snp = tumor2cnv(str(pile), MULTI_BED, str(path), gc=SGC, genmap=SGENMAP)
        if ref is None:
            ref = snp
        pd.testing.assert_frame_equal(snp, ref)
    assert len(ref) >= 20


def test_sites_header_and_record_rules(tmp_path):
    """Only a first line whose Pos has no digit is a TSV header; other lines are data (a bad Pos is an error, never a
    silently skipped site). Multi-allelic sites are dropped as one record or split; identical records count once."""
    cand = pd.DataFrame({"Chr": ["chr7"] * 3, "Pos": [100, 200, 300]})

    def kept(name, text):
        path = tmp_path / name
        path.write_text(text)
        return sorted(sites_at(str(path), cand)[0].itertuples(index=False, name=None))

    assert kept("plus.tsv", "chr7\t+100\nchr7\t200\n") == [("chr7", 100), ("chr7", 200)]  # "+100" is data
    assert kept("hdr.tsv", "#CHROM\tPOS\nchr7\t100\n") == kept("h.tsv", "Chr\tPos\nchr7\t100\n") == [("chr7", 100)]
    with pytest.raises(ValueError, match="'10x' is not a 1-based integer"):
        kept("x.tsv", "chr7\t10x\nchr7\t200\n")
    with pytest.raises(ValueError, match="'Pos' is not a 1-based integer"):  # only the first line can be the header
        kept("hh.tsv", "Chr\tPos\nChr\tPos\nchr7\t100\n")
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n"
    vcf = {"joined": "chr7\t100\t.\tG\tA,T\n", "split": "chr7\t100\t.\tG\tA\nchr7\t100\t.\tG\tT\n",
           "snv+indel": "chr7\t100\t.\tG\tA\nchr7\t100\t.\tG\tGTT\n",
           "snv+star": "chr7\t100\t.\tG\tA\nchr7\t100\t.\tG\t*\n",
           "same twice": "chr7\t100\t.\tG\tA\nchr7\t100\t.\tg\ta\n"}
    for name, rows in vcf.items():
        exp = [("chr7", 100)] if name == "same twice" else []
        assert kept(f"{name}.vcf", head + rows + "chr7\t200\t.\tC\tT\n") == exp + [("chr7", 200)], name
    # a split site at a candidate position is counted in the stats; non-candidate positions are not looked at
    path = tmp_path / "m.vcf"
    path.write_text(head + vcf["split"] + "chr7\t999\t.\tG\tA\nchr7\t999\t.\tG\tC\n")
    assert sites_at(str(path), cand)[1]["multiallelic"] == 1
    # a VCF record that is not kept may have POS 0 (telomere); a kept one may not
    assert kept("tel.vcf", head + "chr7\t0\t.\tN\t<TEL>\nchr7\t200\t.\tC\tT\n") == [("chr7", 200)]
    with pytest.raises(ValueError, match="'0' is not a 1-based integer"):
        kept("zero.vcf", head + "chr7\t0\t.\tC\tT\n")


def test_check_sites(tmp_path):
    """Checked before the pileup pass: missing, not a regular file (read more than once), no record, a BED."""
    for name, text in [("e.tsv", ""), ("h.tsv", "Chr\tPos\n"), ("c.tsv", "#CHROM\tPOS\n"),
                       ("h.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")]:
        (tmp_path / name).write_text(text)
        with pytest.raises(ValueError, match="has no records .empty or header only."):
            check_sites(str(tmp_path / name))
    with pytest.raises(FileNotFoundError, match="does not exist"):
        check_sites(str(tmp_path / "none.vcf"))
    os.mkfifo(tmp_path / "fifo")
    for path in [tmp_path, tmp_path / "fifo"]:
        with pytest.raises(ValueError, match="is not a regular file"):
            check_sites(str(path))
    (tmp_path / "s.bed").write_text("chr7\t99\t100\n")
    with pytest.raises(ValueError, match="2 columns Chr Pos"):
        check_sites(str(tmp_path / "s.bed"))
    (tmp_path / "ok.tsv").write_text("Chr\tPos\nchr7\t100\n")
    assert check_sites(str(tmp_path / "ok.tsv")) == "tsv"


def test_header_only_sites_fail_before_the_pileup_pass(tmp_path, tumor, caplog):
    _, pile, _, _, _ = tumor
    sites = tmp_path / "h.tsv"
    sites.write_text("Chr\tPos\n")
    with caplog.at_level(logging.INFO, logger="mycnv"), pytest.raises(ValueError, match=f"--snp-sites {sites} has no"):
        tumor2cnv(str(pile), MULTI_BED, str(sites), gc=SGC, genmap=SGENMAP)
    assert "Coverage and het SNPs" not in caplog.text


def test_a_split_chromosome_without_sites_is_not_an_error(tmp_path, tumor, caplog):
    """The naming check compares with the BED's chromosomes: a --chrom split with no site on its chromosome (a small
    panel) gives no het SNP and a warning, not an error."""
    _, pile, _, records, _ = tumor
    no_x = tmp_path / "nox.tsv"
    no_x.write_text("".join(f"{c}\t{p}\n" for c, p, *_ in records if c != "chrX"))
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        cov, snp = tumor2cnv(str(pile), MULTI_BED, str(no_x), gc=SGC, genmap=SGENMAP, chrom="chrX")
    assert len(cov) > 0 and snp.empty and "No het SNP" in caplog.text
    _, snp7 = tumor2cnv(str(pile), MULTI_BED, str(no_x), gc=SGC, genmap=SGENMAP, chrom="chr7")
    assert len(snp7) > 0


def test_raw_one_sample_mpileup_route(tmp_path):
    """--raw-pileup with a one-sample raw mpileup (cleanpileup -d added) == the clean route."""
    raw = tmp_path / "t.mpileup"
    with open(f"{SYNTH}/raw.mpileup") as fh:
        raw.write_text("".join("\t".join(line.split("\t")[:6]) + "\n" for line in fh.read().splitlines()))
    clean = tmp_path / "t.cpile"
    with open(clean, "w") as fh:
        subprocess.run(["bash", "-o", "pipefail", "-c", f"{mawk('cleanpileup')} -d < {raw}"], stdout=fh, check=True)
    assert pileup_samples(str(raw), clean=False) == pileup_samples(str(clean)) == 1
    sites = tmp_path / "s.tsv"
    sites.write_text("chr7\t1001\nchr7\t1002\nchr21\t501\n")
    a = tumor2cnv(str(raw), MULTI_BED, str(sites), gc=SGC, genmap=SGENMAP, clean=False, min_depth=4)
    b = tumor2cnv(str(clean), MULTI_BED, str(sites), gc=SGC, genmap=SGENMAP, min_depth=4)
    for x, y in zip(a, b):
        pd.testing.assert_frame_equal(x, y)
    assert a[1]["Pos"].tolist() == [1002, 501]  # chr7:1002 I/i 2 of 6, chr21:501 A/a 2 of 4; 1001 has no alt read


def test_pileup_samples(tmp_path):
    assert pileup_samples(MULTI_PILE) == 2 and pileup_samples(f"{SYNTH}/raw.mpileup", clean=False) == 2
    empty = tmp_path / "e.cpile"
    empty.write_text("")
    with pytest.raises(ValueError, match="is empty"):
        pileup_samples(str(empty))
    with pytest.raises(ValueError, match="needs --raw-pileup"):
        pileup_samples(f"{SYNTH}/raw.mpileup")
    with pytest.raises(ValueError, match="without --raw-pileup"):
        pileup_samples(MULTI_PILE, clean=False)


def test_streamed_pileups(tmp_path, tumor):
    """A pileup on a pipe (/dev/stdin) is read once, by the pass: a pair gives the same files as from a regular file
    (as before the sample-count check), a tumor-only sample too; the sample count and the layout are checked on the
    stream header, each error one line with exit 2 (no traceback)."""
    _, pile, text, _, vcf = tumor
    with open(MULTI_PILE) as fh:
        pair = fh.read()
    with open(f"{SYNTH}/raw.mpileup") as fh:
        raw = fh.read()
    common = ["--bed", MULTI_BED, *TRACKS]

    def tn2cnv_cli(name, pileup, *extra, stdin=None):
        return subprocess.run([*MYCNV, "tn2cnv", "--pileup", str(pileup), *map(str, extra), *common, "--out-cov",
                               str(tmp_path / f"{name}.cov"), "--out-snp", str(tmp_path / f"{name}.snp")],
                              input=stdin, capture_output=True, text=True)

    sites = ["--tumor-only", "--snp-sites", vcf]
    for name, pileup, stdin, extra in [("pair", MULTI_PILE, pair, []), ("tumor", pile, text, sites)]:
        assert tn2cnv_cli(f"{name}_file", pileup, *extra).returncode == 0
        res = tn2cnv_cli(f"{name}_pipe", "/dev/stdin", *extra, stdin=stdin)
        assert res.returncode == 0, res.stderr
        for ext in ["cov", "snp"]:
            with open(tmp_path / f"{name}_file.{ext}", "rb") as a, open(tmp_path / f"{name}_pipe.{ext}", "rb") as b:
                assert a.read() == b.read(), f"{name}.{ext}"
    for stdin, extra, message in [
        (text, [], "/dev/stdin has 1 sample; a pair needs NORMAL TUMOR (2): a one-sample tumor pileup needs"),
        (pair, sites, "--tumor-only needs a one-sample (TUMOR) pileup; /dev/stdin has 2 samples"),
        (raw, [], "/dev/stdin does not start with the `cleanpileup -d` header"),
        ("", [], "the pileup /dev/stdin is empty"),
    ]:
        res = tn2cnv_cli("x", "/dev/stdin", *extra, stdin=stdin)
        assert res.returncode == 2 and message in res.stderr and "Traceback" not in res.stderr, res.stderr
        assert not (tmp_path / "x.cov").exists()


def test_tumor2cnv_rejects_wrong_inputs(tmp_path, tumor):
    _, pile, _, records, vcf = tumor
    kw = dict(gc=SGC, genmap=SGENMAP)
    with pytest.raises(ValueError, match="one-sample .TUMOR. pileup.*run without --tumor-only"):
        tumor2cnv(MULTI_PILE, MULTI_BED, str(vcf), **kw)
    with pytest.raises(ValueError, match="has 1 sample.*--tumor-only --snp-sites"):
        tn2cnv(str(pile), MULTI_BED, **kw)
    with pytest.raises(FileNotFoundError, match="does not exist"):
        tumor2cnv(str(pile), MULTI_BED, str(tmp_path / "missing.vcf.gz"), **kw)
    bed = tmp_path / "sites.bed"
    bed.write_text("".join(f"{c}\t{p - 1}\t{p}\n" for c, p, *_ in records))
    with pytest.raises(ValueError, match="2 columns Chr Pos"):
        tumor2cnv(str(pile), MULTI_BED, str(bed), **kw)
    nochr = tmp_path / "nochr.tsv"
    nochr.write_text("".join(f"{c[3:]}\t{p}\n" for c, p, *_ in records))
    with pytest.raises(ValueError, match="chromosome naming"):
        tumor2cnv(str(pile), MULTI_BED, str(nochr), **kw)
    indels = write_gz(tmp_path / "indels.vcf.gz", synth.sites_vcf([(c, p, "GA", "G") for c, p, *_ in records]))
    with pytest.raises(ValueError, match="no biallelic SNV records"):
        tumor2cnv(str(pile), MULTI_BED, str(indels), **kw)
    bad_pos = tmp_path / "bad.tsv"
    bad_pos.write_text("chr7\t1001\nchr7\t10x\n")
    with pytest.raises(ValueError, match="not a 1-based integer"):
        tumor2cnv(str(pile), MULTI_BED, str(bad_pos), **kw)
    for band in [(0.9, 0.1), (0.1, 1.2), ("a", 1), (0.5,)]:
        with pytest.raises(ValueError, match="tumorOnlyVAF"):
            tumor2cnv(str(pile), MULTI_BED, str(vcf), tumor_vaf=band, **kw)


# ---------- downstream with one sample column: combine, rolling, ASCAT, sex.tsv ---------------------------------
SAMPLE_COLS = re.compile(r"^(Cov|Depth|VAF|absVAF|offVAF|log2ratio|covLLH)([12])(.*)$")


def one_sample(df):
    """The tumor columns of a pair's table as a tumor-only table (*2 -> *1; the normal's *1 dropped)."""
    keep, names = [], {}
    for col in df.columns:
        m = SAMPLE_COLS.match(col)
        if m and m.group(2) == "1":
            continue
        keep.append(col)
        if m:
            names[col] = f"{m.group(1)}1{m.group(3)}"
    return df.loc[:, keep].rename(columns=names)


@pytest.fixture
def pon_dir(tmp_path):
    """main's CNV tables of the synthetic CNV PON as a --ponmatrix folder."""
    pon = tmp_path / "pon"
    pon.mkdir()
    os.symlink(f"{LEGACY}/pon_synth_CNV", pon / "CNV")
    return str(pon)


def synth_sample(tumor_x=1.0):
    files = sorted(os.listdir(f"{SYNTH}/sample"))
    cov = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".cov.gz")])
    snp = concat_tables([f"{SYNTH}/sample/{f}" for f in files if f.endswith(".snp")])
    cov.loc[cov["Chr"] == "chrX", "Cov2"] *= tumor_x
    return cov, snp


def test_one_sample_downstream_equals_the_pairs_tumor(tmp_path, pon_dir, study_config):
    """combine, rolling and the ASCAT tables work per sample column: a tumor-only run gives the pair's tumor values."""
    cov, snp = synth_sample()
    cov_p, snp_p = combine_sample(cov, snp, pon_dir)
    cov_t, snp_t = combine_sample(one_sample(cov), one_sample(snp), pon_dir)
    pd.testing.assert_frame_equal(cov_t, one_sample(cov_p))
    pd.testing.assert_frame_equal(snp_t, one_sample(snp_p))
    roll_p, cnv_p, _ = rolling_step(cov_p, snp_p, study_config)
    roll_t, cnv_t, _ = rolling_step(cov_t, snp_t, study_config)
    assert len(roll_t) > 30
    pd.testing.assert_frame_equal(roll_t, one_sample(roll_p))
    pd.testing.assert_frame_equal(cnv_t, one_sample(cnv_p))
    paired = write_ascat(roll_p, "syn_T-N", str(tmp_path / "p"))
    files = write_ascat(roll_t, "syn_T-N", str(tmp_path / "t"), tumor="T", tumor_only=True)
    assert [os.path.basename(f) for f in files] == ["syn_T-N_T_baf.tsv", "syn_T-N_germline_baf.tsv",
                                                     "syn_T-N_T_logr.tsv", "syn_T-N_germline_logr.tsv"]
    assert sorted(os.listdir(tmp_path / "t")) == sorted(os.path.basename(f) for f in files)
    for path in files:  # byte-identical to the pair's tumor and mock germline tables
        with open(path, "rb") as a, open(tmp_path / "p" / os.path.basename(path), "rb") as b:
            assert a.read() == b.read(), path
    assert len(paired) == 6
    baf, logr = read(files[1]), read(files[3])
    assert list(baf.columns) == ["Unnamed: 0", "chrs", "pos", "syn_T-N"] and len(baf) == len(roll_t)
    assert baf["syn_T-N"].eq(0.5).all() and logr["syn_T-N"].eq(0).all()
    t_baf = read(files[0])
    assert np.allclose(t_baf["syn_T-N"], roll_t["VAF1"]) and set(t_baf["chrs"].astype(str)) == {"7", "21", "X"}


def test_tumor_only_sex_tsv(pon_dir, caplog):
    cfg = default_config()
    for tumor_x, call in [(1.0, "XX"), (0.5, "XY")]:
        cov, snp = synth_sample(tumor_x)
        cov_p, _ = combine_sample(cov, snp, pon_dir)
        cov_t, _ = combine_sample(one_sample(cov), one_sample(snp), pon_dir)
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger="mycnv"):
            t = sex_step(cov_t, "T", None, cfg)
        p = sex_step(cov_p, "T", "N", cfg)
        assert list(t.columns) == sex.SEX_COLUMNS and len(t) == 1
        assert t[["sample", "role", "column", "call", "ascat_gender"]].values.tolist() == [
            ["T", "tumor", "log2ratio1", call, call]]  # the ASCAT gender is the tumor's own call
        cols = ["chrx_ratio", "log2_diff", "n_chrx_bins", "n_auto_bins", "call"]
        pd.testing.assert_frame_equal(t[cols], p.loc[p["role"] == "tumor", cols].reset_index(drop=True))
        assert p["ascat_gender"].eq("XX").all()  # a pair keeps the normal's gender
        assert ("a female tumor that lost one X looks the same" in caplog.text) == (call == "XY")
    with pytest.raises(ValueError, match="a T/N pair.*run without --tumor-only"):
        sex.tumor_sex_table(cov_p, "T")
    with pytest.raises(ValueError, match="no log2ratio2 column.*use --tumor-only"):
        sex.pair_sex_table(cov_t, "T", "N")


def test_write_ascat_tumor_only_rules(tmp_path):
    cnv = pd.DataFrame({"Chr": ["chr7"] * 3, "Pos": [100, 101, 102], "VAF1": 0.4, "log2ratio1_mean": 0.2})
    files = write_ascat(cnv, "P1", str(tmp_path / "a"), tumor_only=True)  # the tumor label defaults to the sample
    assert sorted(os.path.basename(f) for f in files) == ["P1_P1_baf.tsv", "P1_P1_logr.tsv", "P1_germline_baf.tsv",
                                                          "P1_germline_logr.tsv"]
    with open(files[0]) as fh:
        assert fh.readline() == "\tchrs\tpos\tP1\n"
    paired = cnv.assign(VAF2=0.3, log2ratio2_mean=0.1)
    for kw, message in [(dict(normal="N"), "no normal"), (dict(mock_germline=False), "--no-mock-germline"),
                        (dict(tumor="germline"), "clashing")]:
        with pytest.raises(ValueError, match=message):
            write_ascat(cnv, "P1", str(tmp_path / "b"), tumor_only=True, **kw)
    with pytest.raises(ValueError, match="run without --tumor-only"):
        write_ascat(paired, "P1_T-N", str(tmp_path / "b"), tumor_only=True)
    with pytest.raises(ValueError, match="need VAF1/VAF2.*use --tumor-only"):
        write_ascat(cnv, "P1_T-N", str(tmp_path / "b"))
    with pytest.raises(ValueError, match="no SNP passed"):
        write_ascat(cnv.iloc[:0], "P1", str(tmp_path / "b"), tumor_only=True)
    assert not (tmp_path / "b").exists()
    assert len(write_ascat(cnv.iloc[:0], "P1", str(tmp_path / "c"), tumor_only=True, allow_empty=True)) == 4


def test_check_sample_columns():
    pair = pd.DataFrame(columns=["Chr", "Pos", "Cov1", "Cov2"])
    single = pd.DataFrame(columns=["Chr", "Pos", "Depth1", "VAF1"])
    check_sample_columns(False, pair)
    check_sample_columns(True, single)
    with pytest.raises(ValueError, match="use --tumor-only"):
        check_sample_columns(False, single)
    with pytest.raises(ValueError, match="run without --tumor-only"):
        check_sample_columns(True, pair)


def test_tumor_only_plots(tmp_path):
    cov = pd.DataFrame({"Chr": ["chr7"] * 30 + ["chr21"] * 30, "Pos": 10_000 + np.r_[np.arange(30), np.arange(30)] * 10,
                        "ExonPos": np.r_[np.arange(30), np.arange(30)] * 10 + 1, "FullExonPos": np.arange(60) * 10 + 1,
                        "log2ratio1": 0.3, "log2ratio1_mean": 0.3})
    snp = cov.iloc[:5].assign(VAF1=0.6)
    files = cnv_plots(snp, cov.assign(VAF1=np.nan), "S1", "T1", None, str(tmp_path))
    files += raw_plots(cov, snp, "S1", "T1", None, str(tmp_path))
    assert sorted(os.path.basename(f) for f in files) == ["S1_T1.cnv.jpg", "S1_T1.cov.raw.jpg", "S1_T1.snp.cnv.jpg",
                                                          "S1_T1.snp.raw.jpg"]
    assert all(os.path.getsize(f) > 0 for f in files)


def test_config_band():
    assert default_config()["hetSNP"]["tumorOnlyVAF"] == [0.1, 0.9]
    assert load_config(overrides={"hetSNP": {"tumorOnlyVAF": [0.2, 0.8]}})["hetSNP"]["tumorOnlyVAF"] == [0.2, 0.8]
    assert "tumorOnlyVAF:\n    - 0.1\n    - 0.9" in run("config").stdout


# ---------- the CLI -----------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def multi_pon(tmp_path_factory):
    pon = tmp_path_factory.mktemp("multi_to") / "PON"
    run("pon-coverage", "--matrix", *MULTI_MATRIX, "--bed", MULTI_BED, *TRACKS, "-o", pon)
    run("pon-combine", "--ponmatrix", pon, "-o", pon, "--chrom", "chr7", "--chrom", "chr21", "--chrom", "chrX")
    return pon


def test_sample_tumor_only_end_to_end_and_split_path(tmp_path, tumor, multi_pon):
    """`sample --tumor-only` writes one-sample tables, a one-row sex.tsv, the tumor + mock germline tables and the
    tumor's plots; tn2cnv --tumor-only per chromosome -> combine -> rolling -> ascat-input gives the same files."""
    _, pile, _, _, vcf = tumor
    out = tmp_path / "s"
    res = run("sample", "--tumor-only", "--snp-sites", vcf, "--pileup", pile, "--bed", MULTI_BED, *TRACKS,
              "--ponmatrix", multi_pon, "--sample", "M_T", "-o", out, "--plots")  # SNPs pass: no --allow-empty
    raw = read(out / "raw" / "M_T.cov.gz")
    assert {"Cov1", "log2ratio1", "PONcov_mean"} <= set(raw.columns) and "Cov2" not in raw.columns
    assert raw["log2ratio1"].notna().mean() > 0.9
    snp = read(out / "raw" / "M_T.snp")
    assert len(snp) >= 20 and "VAF2" not in snp.columns and snp["VAF1"].between(*BAND).all()
    assert {"log2ratio1_mean", "VAF1"} <= set(read(out / "data" / "M_T.cnv.gz").columns)
    assert sorted(os.listdir(out / "pre")) == ["M_T_M_T_baf.tsv", "M_T_M_T_logr.tsv", "M_T_germline_baf.tsv",
                                               "M_T_germline_logr.tsv"]  # the tumor label defaults to --sample
    baf = read(out / "pre" / "M_T_M_T_baf.tsv")
    assert len(baf) > 0 and list(baf.columns) == ["Unnamed: 0", "chrs", "pos", "M_T"]
    assert read(out / "pre" / "M_T_germline_baf.tsv")["M_T"].eq(0.5).all()
    assert read(out / "pre" / "M_T_germline_logr.tsv")["M_T"].eq(0).all()
    s = read(out / "sex.tsv")
    assert s[["sample", "role", "column"]].values.tolist() == [["M_T", "tumor", "log2ratio1"]]
    assert s["call"].tolist() == ["unknown"] and s["ascat_gender"].tolist() == ["XX"]  # few chrX bins
    assert "Sex call of the tumor M_T is unknown" in res.stderr
    plots = sorted(os.listdir(out / "plots"))
    assert plots == ["M_T_M_T.cnv.jpg", "M_T_M_T.cov.raw.jpg", "M_T_M_T.snp.cnv.jpg", "M_T_M_T.snp.raw.jpg"]
    # --tumor-name
    run("sample", "--tumor-only", "--snp-sites", vcf, "--pileup", pile, "--bed", MULTI_BED, *TRACKS, "--ponmatrix",
        multi_pon, "--sample", "M_T", "--tumor-name", "TUM", "-o", tmp_path / "named", "--allow-empty")
    assert sorted(os.listdir(tmp_path / "named" / "pre"))[0] == "M_T_TUM_baf.tsv"
    assert read(tmp_path / "named" / "sex.tsv")["sample"].tolist() == ["TUM"]
    # split path
    covs, snps = [], []
    for chrom in ["chr7", "chr21", "chrX"]:
        covs.append(tmp_path / f"{chrom}.cov.gz")
        snps.append(tmp_path / f"{chrom}.snp")
        run("tn2cnv", "--tumor-only", "--snp-sites", vcf, "--pileup", pile, "--bed", MULTI_BED, *TRACKS,
            "--chrom", chrom, "--out-cov", covs[-1], "--out-snp", snps[-1])
    run("combine", "--tumor-only", "--cov", *covs, "--snp", *snps, "--ponmatrix", multi_pon, "--out-cov",
        tmp_path / "raw.cov.gz", "--out-snp", tmp_path / "raw.snp", "--sample", "M_T")
    run("sex", "--tumor-only", "--cov", tmp_path / "raw.cov.gz", "-o", tmp_path / "sex2.tsv", "--sample", "M_T")
    run("rolling", "--tumor-only", "--cov", tmp_path / "raw.cov.gz", "--snp", tmp_path / "raw.snp", "--out-snp",
        tmp_path / "d.cnv.snp.gz", "--out-cnv", tmp_path / "d.cnv.gz", "--plots-dir", tmp_path / "plots",
        "--sample", "M_T")
    run("ascat-input", "--tumor-only", "--cnv", tmp_path / "d.cnv.snp.gz", "--sample", "M_T", "-o", tmp_path / "pre")
    for split, ref in [("raw.cov.gz", "raw/M_T.cov.gz"), ("raw.snp", "raw/M_T.snp"), ("d.cnv.gz", "data/M_T.cnv.gz"),
                       ("d.cnv.snp.gz", "data/M_T.cnv.snp.gz"), ("sex.tsv", "sex.tsv"), ("sex2.tsv", "sex.tsv")]:
        pd.testing.assert_frame_equal(read(tmp_path / split), read(out / ref))
    for name in os.listdir(out / "pre"):
        with open(tmp_path / "pre" / name) as a, open(out / "pre" / name) as b:
            assert a.read() == b.read()
    assert sorted(os.listdir(tmp_path / "plots")) == ["M_T_M_T.cnv.jpg", "M_T_M_T.snp.cnv.jpg"]


def test_step_by_step_tumor_only_with_snps_equals_the_pairs_tumor(tmp_path, pon_dir):
    """combine/sex/rolling/ascat-input --tumor-only on one-sample tables (the synthetic sample's tumor columns): real
    BAF rows, the same tumor and mock germline tables as the pair's chain."""
    files = sorted(os.listdir(f"{SYNTH}/sample"))
    cov_in, snp_in = [], []
    for f in files:
        path = tmp_path / f"t_{f}"
        one_sample(read(f"{SYNTH}/sample/{f}")).to_csv(path, sep="\t", index=False)
        (cov_in if f.endswith(".cov.gz") else snp_in).append(path)
    pair = {"cov": [f"{SYNTH}/sample/{f}" for f in files if f.endswith(".cov.gz")],
            "snp": [f"{SYNTH}/sample/{f}" for f in files if f.endswith(".snp")]}
    for tag, cov, snp, extra in [("t", cov_in, snp_in, ["--tumor-only"]), ("p", pair["cov"], pair["snp"], [])]:
        d = tmp_path / tag
        run("combine", *extra, "--cov", *cov, "--snp", *snp, "--ponmatrix", pon_dir, "--out-cov", d / "raw.cov.gz",
            "--out-snp", d / "raw.snp", "--tumor-name", "T", *([] if extra else ["--normal-name", "N"]))
        run("rolling", *extra, "--cov", d / "raw.cov.gz", "--snp", d / "raw.snp", "--out-snp", d / "x.cnv.snp.gz",
            "--out-cnv", d / "x.cnv.gz")
        run("ascat-input", *extra, "--cnv", d / "x.cnv.snp.gz", "--sample", "syn_T-N", "-o", d / "pre",
            "--tumor-name", "T", *([] if extra else ["--normal-name", "N"]))
    t_files = sorted(os.listdir(tmp_path / "t" / "pre"))
    assert t_files == ["syn_T-N_T_baf.tsv", "syn_T-N_T_logr.tsv", "syn_T-N_germline_baf.tsv",
                       "syn_T-N_germline_logr.tsv"]
    for name in t_files:
        with open(tmp_path / "t" / "pre" / name) as a, open(tmp_path / "p" / "pre" / name) as b:
            assert a.read() == b.read()
    assert len(read(tmp_path / "t" / "pre" / "syn_T-N_T_baf.tsv")) > 30
    s = read(tmp_path / "t" / "sex.tsv")  # combine's default: next to --out-cov
    assert s[["sample", "role", "call", "ascat_gender"]].values.tolist() == [["T", "tumor", "XX", "XX"]]
    # tumor-only and paired tables are not interchangeable
    res = run("combine", "--cov", *cov_in, "--snp", *snp_in, "--ponmatrix", pon_dir, "--out-cov", tmp_path / "x" / "c",
              "--out-snp", tmp_path / "x" / "s", check=False)
    assert res.returncode == 2 and "use --tumor-only" in res.stderr and not (tmp_path / "x").exists()
    res = run("combine", "--tumor-only", *["--cov", *pair["cov"]], "--snp", *pair["snp"], "--ponmatrix", pon_dir,
              "--out-cov", tmp_path / "x" / "c", "--out-snp", tmp_path / "x" / "s", check=False)
    assert res.returncode == 2 and "run without --tumor-only" in res.stderr and not (tmp_path / "x").exists()
    res = run("ascat-input", "--tumor-only", "--cnv", tmp_path / "p" / "x.cnv.snp.gz", "--sample", "S", "-o",
              tmp_path / "x", check=False)
    assert res.returncode == 2 and "run without --tumor-only" in res.stderr
    res = run("rolling", "--cov", tmp_path / "t" / "raw.cov.gz", "--snp", tmp_path / "t" / "raw.snp", "--out-snp",
              tmp_path / "x" / "a", "--out-cnv", tmp_path / "x" / "b", check=False)
    assert res.returncode == 2 and "use --tumor-only" in res.stderr
    res = run("sex", "--tumor-only", "--cov", tmp_path / "p" / "raw.cov.gz", "-o", tmp_path / "x" / "s.tsv", check=False)
    assert res.returncode == 2 and "run without --tumor-only" in res.stderr and not (tmp_path / "x").exists()


@pytest.mark.parametrize("argv, message", [
    (["sample", "--tumor-only"], "--tumor-only needs --snp-sites"),
    (["sample", "--snp-sites", "s.vcf"], "--snp-sites goes with --tumor-only"),
    (["sample", "--tumor-only", "--snp-sites", "s.vcf", "--normal-name", "N"], "--normal-name does not go with"),
    (["sample", "--tumor-only", "--snp-sites", "s.vcf", "--no-mock-germline"], "--no-mock-germline does not go with"),
    (["sample", "--tumor-only", "--snp-sites", "s.vcf", "--normal-vaf", "0.2", "0.8"], "--normal-vaf does not apply"),
    (["sample", "--tumor-vaf", "0.2", "0.8"], "--tumor-vaf goes with --tumor-only"),
    (["sample", "--tumor-name", "T"], "give both --tumor-name and --normal-name"),  # pairs: unchanged
    (["tn2cnv", "--tumor-only"], "--tumor-only needs --snp-sites"),
    (["tn2cnv", "--snp-sites", "s.vcf"], "--snp-sites goes with --tumor-only"),
    (["ascat-input", "--tumor-only", "--no-mock-germline"], "--no-mock-germline does not go with"),
    (["sex", "--tumor-only", "--normal-name", "N"], "--normal-name does not go with"),
])
def test_rejected_tumor_only_arguments(tmp_path, argv, message):
    base = {
        "sample": ["--pileup", "p", "--ponmatrix", "d", "--sample", "s", "-o", tmp_path / "o", "--bed", "b", *TRACKS],
        "tn2cnv": ["--pileup", "p", "--bed", "b", *TRACKS, "--out-cov", tmp_path / "c", "--out-snp", tmp_path / "s"],
        "ascat-input": ["--cnv", "x", "--sample", "s", "-o", tmp_path / "o"],
        "sex": ["--cov", "x", "-o", tmp_path / "s.tsv"],
    }[argv[0]]
    res = run(argv[0], *base, *argv[1:], check=False)
    assert res.returncode == 2 and "error:" in res.stderr and message in res.stderr, res.stderr


def test_sample_rejects_the_other_pileup_layout(tmp_path, tumor, multi_pon):
    _, pile, _, _, vcf = tumor
    common = ["--bed", MULTI_BED, *TRACKS, "--ponmatrix", multi_pon, "--sample", "M_T-N", "--allow-empty"]
    res = run("sample", "--tumor-only", "--snp-sites", vcf, "--pileup", MULTI_PILE, "-o", tmp_path / "a", *common,
              check=False)
    assert res.returncode == 2 and "one-sample (TUMOR) pileup" in res.stderr and "run without --tumor-only" in res.stderr
    res = run("sample", "--pileup", pile, "-o", tmp_path / "b", *common, check=False)
    assert res.returncode == 2 and "needs --tumor-only --snp-sites" in res.stderr
    cfg = tmp_path / "c.yaml"
    cfg.write_text("CNV:\n  hetSNP:\n    tumorOnlyVAF: [0.9, 0.1]\n")
    res = run("sample", "--tumor-only", "--snp-sites", vcf, "--pileup", pile, "-o", tmp_path / "c", *common,
              "--config", cfg, check=False)
    assert res.returncode == 2 and "tumorOnlyVAF needs 0 <= min <= max <= 1" in res.stderr
    res = run("sample", "--tumor-only", "--snp-sites", tmp_path / "none.vcf.gz", "--pileup", pile, "-o",
              tmp_path / "d", *common, check=False)
    assert res.returncode == 2 and "does not exist" in res.stderr
    (tmp_path / "h.tsv").write_text("Chr\tPos\n")
    res = run("sample", "--tumor-only", "--snp-sites", tmp_path / "h.tsv", "--pileup", pile, "-o", tmp_path / "e",
              *common, check=False)
    assert res.returncode == 2 and f"--snp-sites {tmp_path / 'h.tsv'} has no records" in res.stderr
    assert not any((tmp_path / x).exists() for x in "abcde")
    # --tumor-vaf reaches the band
    run("tn2cnv", "--tumor-only", "--snp-sites", vcf, "--pileup", pile, "--bed", MULTI_BED, *TRACKS, "--out-cov",
        tmp_path / "n.cov", "--out-snp", tmp_path / "n.snp", "--tumor-vaf", 0.4, 0.6)
    n = read(tmp_path / "n.snp")
    assert len(n) > 0 and n["VAF1"].between(0.4, 0.6).all()
