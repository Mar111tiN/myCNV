"""Bin coordinates shared by samples and PONs (fixed mode), BED checks, and the sample/PON compatibility guards."""

import gzip
import os
import subprocess
import sys
from io import StringIO

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, GC, GENMAP, SYNTH, read
from mycnv.bed import bin_exon_pos, read_bed, region_digest
from mycnv.raw import pon2cnv, tn2cnv
from mycnv.tools import mawk

BED = f"{DATA}/haev7_covered_chr7.bed"
PADDED = f"{DATA}/haev7_padded_chr7.bed"
MYCNV = [sys.executable, "-m", "mycnv.cli"]


def sh(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)


def cli(*args, check=True):
    res = subprocess.run([*MYCNV, *map(str, args)], capture_output=True, text=True)
    if check and res.returncode:
        raise AssertionError(f"mycnv {' '.join(map(str, args))} failed:\n{res.stderr[-3000:]}")
    return res


def dense_first_base_exonpos(bed: str, chrom: str = "chr7") -> pd.Series:
    """ExonPos of the first filterBed-kept base of every 10 bp bin, from the real tool on a dense position stream."""
    end = max(read_bed(bed)[chrom][1]) + 20
    stream = "Chr\tStart\n" + "".join(f"{chrom}\t{p}\n" for p in range(1, end))
    res = subprocess.run(["bash", "-o", "pipefail", "-c", f"{mawk('filterBed')} {bed} -x"], input=stream,
                         capture_output=True, text=True, check=True)
    df = pd.read_csv(StringIO(res.stdout), sep="\t")
    df["Pos"] = df["Start"] // 10 * 10
    return df.drop_duplicates("Pos").set_index("Pos")["ExonPos"]


def thinned_pileup(path, drop_every=7):
    """chr7.cpile without the first two bases of every BED region and without every `drop_every`-th row, i.e. a
    stream with other positions than the PONmatrix (as a `samtools mpileup -l` pileup vs a padded-BED matrix)."""
    pile = read(f"{SYNTH}/chr7.cpile")
    starts = set()
    for s in read_bed(BED)["chr7"][0]:
        starts.update({int(s), int(s) + 1})
    keep = ~pile["Start"].isin(starts) & (np.arange(len(pile)) % drop_every != 3)
    pile.loc[keep].to_csv(path, sep="\t", index=False)
    return path


# ---------- BED handling ----------------------------------------------------------------------------------------
@pytest.mark.parametrize("lines, message", [
    (["chr1\t2009\t2030", "chr1\t1000\t1001"], "listed after"),
    (["chr1\t100\t105", "chr2\t50\t80", "chr1\t2009\t2030"], "comes back"),
])
def test_unsorted_or_interleaved_bed_is_an_error(tmp_path, lines, message):
    bed = tmp_path / "bad.bed"
    bed.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match=message):
        read_bed(str(bed))
    stream = tmp_path / "stream"
    stream.write_text("Chr\tStart\nchr1\t101\nchr1\t1000\nchr1\t2010\n")
    res = sh(f"{mawk('filterBed')} {bed} -x < {stream}")
    assert res.returncode != 0 and message.split()[0] in res.stderr


def test_filterbed_chromosome_list():
    res = sh(f"{mawk('filterBed')} {SYNTH}/multi.bed -x -c 21,chrX < {SYNTH}/multi.cpile")
    df = pd.read_csv(StringIO(res.stdout), sep="\t")
    assert set(df["Chr"]) == {"chr21", "chrX"}


def test_bin_exonpos_is_the_first_bed_base_of_the_bin():
    expected = dense_first_base_exonpos(BED)
    cov, _ = tn2cnv(f"{SYNTH}/chr7.cpile", BED, gc=GC, genmap=GENMAP)
    assert (cov["Pos"] % 10 == 0).all()
    assert np.array_equal(cov["ExonPos"].to_numpy(), expected.loc[cov["Pos"]].to_numpy())
    # GC is merged on the bin start (own label)
    gc = pd.read_csv(f"{GC}/chr7.gc100-10.gz", sep="\t").set_index("Start")["GC/AT"]
    assert np.allclose(cov["GCratio"], gc.loc[cov["Pos"]].to_numpy())


def test_bin_exonpos_with_overlapping_and_bookended_regions(tmp_path):
    bed = tmp_path / "edge.bed"
    bed.write_text("chr7\t100\t200\nchr7\t150\t180\nchr7\t203\t205\nchr7\t205\t230\nchr7\t236\t238\nchr7\t300\t301\n")
    expected = dense_first_base_exonpos(str(bed))
    bins = pd.DataFrame({"Chr": "chr7", "Pos": expected.index})
    assert np.array_equal(bin_exon_pos(bins, read_bed(str(bed))), expected.to_numpy())


OVERLAP_BED = ("chr7\t100\t200\nchr7\t150\t180\nchr7\t190\t260\nchr7\t260\t270\nchr7\t271\t280\n"
               "chr7\t400\t410\nchr8\t10\t50\nchr8\t20\t30\n")


def _filterbed_dense(bed, legacy=False, end=500):
    stream = "Chr\tStart\n" + "".join(f"{c}\t{p}\n" for c in ("chr7", "chr8") for p in range(1, end))
    res = subprocess.run(["bash", "-o", "pipefail", "-c", f"{mawk('filterBed')} {bed} -x{' -L' if legacy else ''}"],
                         input=stream, capture_output=True, text=True, check=True)
    return pd.read_csv(StringIO(res.stdout), sep="\t")


def test_overlapping_bed_regions_never_double_count_a_base(tmp_path):
    """Vendor padded BEDs overlap (SS_HAEv7 padded: 29,488 overlaps): in the fixed coordinates every kept base gets
    the next ExonPos (1, 2, 3, ... per chromosome), whether the BED is merged or not."""
    bed = tmp_path / "overlap.bed"
    bed.write_text(OVERLAP_BED)
    out = _filterbed_dense(bed)
    for chrom, df in out.groupby("Chr"):
        assert df["Start"].is_unique
        assert df["ExonPos"].tolist() == list(range(1, len(df) + 1)), chrom
    # the kept positions are the union of the regions' start..end (filterBed convention, both inclusive)
    exp = {("chr7", p) for s, e in [(100, 280), (400, 410)] for p in range(s, e + 1)}
    exp |= {("chr8", p) for p in range(10, 51)}
    assert set(zip(out["Chr"], out["Start"])) == exp


def test_overlapping_bed_equals_its_merged_copy(tmp_path):
    raw = tmp_path / "overlap.bed"
    raw.write_text(OVERLAP_BED)
    merged = read_bed(str(raw))
    assert {c: (s.tolist(), e.tolist()) for c, (s, e) in merged.items()} == {
        "chr7": ([100, 400], [280, 410]), "chr8": ([10], [50])}
    mfile = tmp_path / "merged.bed"
    mfile.write_text("".join(f"{c}\t{a}\t{b}\n" for c, (s, e) in merged.items() for a, b in zip(s, e)))
    pd.testing.assert_frame_equal(_filterbed_dense(raw), _filterbed_dense(mfile))
    assert region_digest(read_bed(str(raw)), "chr7") == region_digest(read_bed(str(mfile)), "chr7")
    # BED-derived bin ExonPos (fixed) = ExonPos of the first kept base of the bin, on the raw overlapping BED
    dense = _filterbed_dense(raw)
    dense = dense[dense["Chr"] == "chr7"].assign(Pos=lambda d: d["Start"] // 10 * 10).drop_duplicates("Pos")
    bins = pd.DataFrame({"Chr": "chr7", "Pos": dense["Pos"].to_numpy()})
    assert np.array_equal(bin_exon_pos(bins, read_bed(str(raw))), dense["ExonPos"].to_numpy())


def test_legacy_coords_keep_the_unmerged_regions(tmp_path):
    """-L reproduces main: overlapping regions stay separate, so ExonPos (BEDSUM - BEDSTART + pos) counts them twice."""
    bed = tmp_path / "overlap.bed"
    bed.write_text(OVERLAP_BED)
    out = _filterbed_dense(bed, legacy=True)
    s, e = read_bed(str(bed), merge=False)["chr7"]
    assert s.tolist() == [100, 150, 190, 260, 271, 400]  # unmerged
    bedsum = np.cumsum(e - s + 1)
    row = out[(out["Chr"] == "chr7") & (out["Start"] == 400)].iloc[0]
    assert row["ExonPos"] == bedsum[5] - s[5] + 400  # region 6 of the unmerged list
    fixed = _filterbed_dense(bed)
    assert fixed.loc[(fixed["Chr"] == "chr7") & (fixed["Start"] == 400), "ExonPos"].item() < row["ExonPos"]


def test_tn2cnv_and_pon2cnv_on_an_overlapping_bed_equal_its_merged_copy(tmp_path):
    raw = PADDED  # the SS_HAEv7 padded slice has overlapping regions, like the full kit BED
    s, e = read_bed(raw, merge=False)["chr7"]
    assert (s[1:] <= e[:-1]).any(), "the fixture must contain overlapping regions"
    merged = tmp_path / "merged.bed"
    m = read_bed(raw)
    merged.write_text("".join(f"chr7\t{a}\t{b}\n" for a, b in zip(*m["chr7"])))
    a_cov, a_snp = tn2cnv(f"{SYNTH}/chr7.cpile", raw, gc=GC, genmap=GENMAP)
    b_cov, b_snp = tn2cnv(f"{SYNTH}/chr7.cpile", str(merged), gc=GC, genmap=GENMAP)
    pd.testing.assert_frame_equal(a_cov, b_cov)
    pd.testing.assert_frame_equal(a_snp, b_snp)
    assert a_cov.groupby("Chr")["ExonPos"].apply(lambda x: x.is_monotonic_increasing).all()
    p_cov, _ = pon2cnv(f"{SYNTH}/matrix/chr7.pon.gz", raw, gc=GC, chrom="chr7")
    q_cov, _ = pon2cnv(f"{SYNTH}/matrix/chr7.pon.gz", str(merged), gc=GC, chrom="chr7")
    pd.testing.assert_frame_equal(p_cov, q_cov)


def test_pile2cnv_fixed_bin_values_from_the_stream(tmp_path):
    """Pos = 10 * floor(first base / 10), Cov = sum(depth) / 10 per bin and sample (raw tool output)."""
    snp = tmp_path / "s.snp"
    stream = sh(f"{mawk('filterBed')} {BED} -x < {SYNTH}/chr7.cpile").stdout
    cov = pd.read_csv(StringIO(sh(f"{mawk('filterBed')} {BED} -x < {SYNTH}/chr7.cpile | "
                                  f"{mawk('pile2CNV')} -x -o {snp}").stdout), sep="\t")
    st = pd.read_csv(StringIO(stream), sep="\t")
    st["Pos"] = st["Start"] // 10 * 10
    g = st.groupby("Pos", sort=False)
    exp = pd.DataFrame({"Pos": g["Start"].first() // 10 * 10, "ExonPos": g["ExonPos"].first(),
                        "Cov1": g["Depth1"].sum() / 10, "Cov2": g["Depth2"].sum() / 10}).reset_index(drop=True)
    assert np.array_equal(cov["Pos"], exp["Pos"]) and np.array_equal(cov["ExonPos"], exp["ExonPos"])
    assert np.allclose(cov[["Cov1", "Cov2"]], exp[["Cov1", "Cov2"]])


# ---------- sample vs PON ---------------------------------------------------------------------------------------
def test_sample_and_pon_bins_agree_when_streams_hold_other_positions(tmp_path):
    """Review finding: the first-streamed-base ExonPos lost ~1 bin per region start in the sample/PON join."""
    pile = thinned_pileup(tmp_path / "thin.cpile")
    cov, _ = tn2cnv(str(pile), BED, gc=GC, genmap=GENMAP)
    pon, _ = pon2cnv(f"{SYNTH}/matrix/chr7.pon.gz", BED, gc=GC, genmap=GENMAP)
    m = cov.merge(pon, on=["Chr", "Pos"], suffixes=("", "_pon"))
    assert len(m) == len(cov)
    assert (m["ExonPos"] == m["ExonPos_pon"]).all()


@pytest.fixture(scope="module")
def pons(tmp_path_factory):
    """CNV PONs from the synthetic chr7 PONmatrix: fixed and legacy coordinates, covered and padded BED."""
    root = tmp_path_factory.mktemp("pons")
    out = {}
    for name, bed, extra in [("fixed", BED, []), ("legacy", BED, ["--legacy-coords"]), ("padded", PADDED, [])]:
        cli("pon-coverage", "--ponmatrix", SYNTH, "--bed", bed, "--gc", GC, "--genmap", GENMAP,
            "-o", root / name, *extra)
        cli("pon-combine", "--ponmatrix", root / name, "-o", root / name, "--chrom", "chr7")
        out[name] = root / name
    return out


def test_thinned_sample_keeps_every_bin_through_the_cli(tmp_path, pons):
    pile = thinned_pileup(tmp_path / "thin.cpile")
    cli("sample", "--pileup", pile, "--bed", BED, "--gc", GC, "--genmap", GENMAP, "--ponmatrix", pons["fixed"],
        "--sample", "S1_T-N", "-o", tmp_path / "out", "--allow-empty")
    raw = read(tmp_path / "out" / "raw" / "S1_T-N.cov.gz")
    cov, _ = tn2cnv(str(pile), BED, gc=GC, genmap=GENMAP)
    pon = read(pons["fixed"] / "CNV" / "pon.cov.filterGC.gz")
    shared = cov[["Pos"]].merge(pon[["Pos"]])
    assert len(raw) == len(shared) == len(cov)


def test_pon_records_its_coordinates(pons):
    import json

    meta = json.loads((pons["legacy"] / "CNV" / "pon.json").read_text())
    assert meta["coords"] == "legacy" and meta["chroms"] == ["chr7"] and meta["bin_size"] == 10
    assert json.loads((pons["fixed"] / "cov" / "chr7.cov.json").read_text())["coords"] == "fixed"


@pytest.mark.parametrize("pon, extra, message", [
    ("legacy", [], "add --legacy-coords"),
    ("fixed", ["--legacy-coords"], "drop --legacy-coords"),
    ("padded", [], "differs from the BED"),
])
def test_mixed_coordinates_or_beds_fail(tmp_path, pons, pon, extra, message):
    res = cli("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED, "--gc", GC, "--genmap", GENMAP,
              "--ponmatrix", pons[pon], "--sample", "S1_T-N", "-o", tmp_path / "out", *extra, check=False)
    assert res.returncode != 0 and message in res.stderr
    assert not (tmp_path / "out" / "raw").exists()


@pytest.mark.parametrize("pon, extra, message", [
    ("legacy", [], "duplicate bin labels|another ExonPos|PON bin"),
    ("padded", [], "another ExonPos"),
])
def test_mismatch_without_pon_metadata_is_caught_by_the_data(tmp_path, pons, pon, extra, message):
    """PONs built by the legacy tools carry no pon.json: the bin keys themselves must reveal the mismatch."""
    import re
    import shutil

    bare = tmp_path / "bare"
    shutil.copytree(pons[pon] / "CNV", bare / "CNV")
    (bare / "CNV" / "pon.json").unlink()
    res = cli("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED, "--gc", GC, "--genmap", GENMAP,
              "--ponmatrix", bare, "--sample", "S1_T-N", "-o", tmp_path / "out", *extra, check=False)
    assert res.returncode != 0 and re.search(message, res.stderr)


def test_legacy_sample_against_fixed_pon_without_metadata_fails(tmp_path, pons):
    import shutil

    bare = tmp_path / "bare"
    shutil.copytree(pons["fixed"] / "CNV", bare / "CNV")
    (bare / "CNV" / "pon.json").unlink()
    res = cli("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED, "--gc", GC, "--genmap", GENMAP,
              "--ponmatrix", bare, "--sample", "S1_T-N", "-o", tmp_path / "out", "--legacy-coords", check=False)
    assert res.returncode != 0 and "have a PON bin" in res.stderr


def test_too_few_pon_bins_fail(tmp_path, pons):
    import shutil

    small = tmp_path / "small"
    shutil.copytree(pons["fixed"] / "CNV", small / "CNV")
    path = small / "CNV" / "pon.cov.filterGC.gz"
    df = read(path)
    df.iloc[: len(df) // 2].to_csv(path, sep="\t", index=False, compression="gzip")
    res = cli("sample", "--pileup", f"{SYNTH}/chr7.cpile", "--bed", BED, "--gc", GC, "--genmap", GENMAP,
              "--ponmatrix", small, "--sample", "S1_T-N", "-o", tmp_path / "out", check=False)
    assert res.returncode != 0 and "have a PON bin" in res.stderr


def test_pon2cnv_fixed_bin_values_from_the_stream(tmp_path):
    """PON2CNV: Pos = 10 * floor(first base / 10), Cov = per-sample (+ and - strand) depth sum / 10; rows without
    coverage in every normal are skipped."""
    snp = tmp_path / "p.snp"
    filt = f"gunzip -c {SYNTH}/matrix/chr7.pon.gz | {mawk('filterBed')} {BED} -x"
    st = pd.read_csv(StringIO(sh(filt).stdout), sep="\t")
    cov = pd.read_csv(StringIO(sh(f"{filt} | {mawk('PON2CNV')} -x -o {snp}").stdout), sep="\t")
    depth = st["Depth"].str.split("=", expand=True)
    d = depth[0].str.split("|", expand=True).astype(int) + depth[1].str.split("|", expand=True).astype(int)
    keep = d.sum(axis=1) > 0
    st, d = st.loc[keep].reset_index(drop=True), d.loc[keep].reset_index(drop=True)
    d.columns = [f"Cov{i + 1}" for i in d.columns]
    d["Pos"] = st["Start"] // 10 * 10
    g = d.groupby("Pos", sort=False)
    exp = (g.sum() / 10).reset_index()
    assert np.array_equal(cov["Pos"], exp["Pos"])
    assert np.array_equal(cov["ExonPos"], st.groupby(d["Pos"], sort=False)["ExonPos"].first().to_numpy())
    cols = [c for c in cov.columns if c.startswith("Cov")]
    assert np.allclose(cov[cols], exp[cols])
