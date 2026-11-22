"""mawk tools: legacy parity (-L) and the fixed semantics (bins, ExonPos, chromosome handling, bash)."""

import os
import subprocess
from io import StringIO

import numpy as np
import pandas as pd

from conftest import DATA, LEGACY, RESTRICTED, RLEGACY, SYNTH, needs_restricted, read
from mycnv.tools import mawk


def run(cmd: str) -> str:
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


def pile2cnv(tmp_path, legacy: bool):
    snp = tmp_path / ("leg.snp" if legacy else "fix.snp")
    flag = " -L" if legacy else ""
    cov = run(f"{mawk('filterBed')} {SYNTH}/multi.bed -x{flag} < {SYNTH}/multi.cpile | "
              f"{mawk('pile2CNV')} -x{flag} -o {snp} -v 0.2 -V 0.8 -d 15")
    return pd.read_csv(StringIO(cov), sep="\t", dtype=str), pd.read_csv(snp, sep="\t", dtype=str)


def test_pile2cnv_legacy_mode_equals_legacy_per_chromosome(tmp_path):
    cov, snp = pile2cnv(tmp_path, legacy=True)
    pd.testing.assert_frame_equal(cov, read(f"{LEGACY}/pile2cnv_multi.cov", dtype=str))
    pd.testing.assert_frame_equal(snp, read(f"{LEGACY}/pile2cnv_multi.snp", dtype=str))


def test_pile2cnv_fixed_bins_and_exonpos(tmp_path):
    cov, snp = pile2cnv(tmp_path, legacy=False)
    leg_cov, leg_snp = pile2cnv(tmp_path, legacy=True)
    cov = cov.astype({"Pos": int, "ExonPos": int})
    # same bins and values, only the labels change
    assert len(cov) == len(leg_cov)
    assert (cov[["Chr", "Cov1", "Cov2"]].to_numpy() == leg_cov[["Chr", "Cov1", "Cov2"]].to_numpy()).all()
    # every bin is labelled with its own start and labels are unique per chromosome
    assert (cov["Pos"] % 10 == 0).all()
    assert not cov.duplicated(["Chr", "Pos"]).any()
    for _, c in cov.groupby("Chr"):
        assert c["Pos"].is_monotonic_increasing and c["ExonPos"].is_monotonic_increasing
    # het-SNP rows are unchanged apart from ExonPos, which is contiguous and increasing per chromosome
    pd.testing.assert_frame_equal(snp.drop(columns="ExonPos"), leg_snp.drop(columns="ExonPos"))
    for _, s in snp.astype({"ExonPos": int}).groupby("Chr"):
        assert (np.diff(s["ExonPos"].to_numpy()) > 0).all()


def test_filterbed_exonpos_is_contiguous():
    out = run(f"{mawk('filterBed')} {SYNTH}/multi.bed -x < {SYNTH}/multi.cpile")
    df = pd.read_csv(StringIO(out), sep="\t")
    for _, c in df.groupby("Chr", sort=False):
        assert c["ExonPos"].iloc[0] == 1
        assert (np.diff(c["ExonPos"].to_numpy()) == 1).all()


def test_filterbed_chromosome_order_and_no_early_exit(tmp_path):
    # stream with chrX first and chr7 last; the BED lists chr7 first and chrX last
    pile = read(f"{SYNTH}/multi.cpile")
    order = {"chrX": 0, "chr21": 1, "chr7": 2}
    shuffled = pile.sort_values("Chr", key=lambda s: s.map(order), kind="stable")
    path = tmp_path / "shuffled.cpile"
    shuffled.to_csv(path, sep="\t", index=False)
    a = pd.read_csv(StringIO(run(f"{mawk('filterBed')} {SYNTH}/multi.bed -x < {SYNTH}/multi.cpile")), sep="\t")
    b = pd.read_csv(StringIO(run(f"{mawk('filterBed')} {SYNTH}/multi.bed -x < {path}")), sep="\t")
    key = ["Chr", "Start"]
    pd.testing.assert_frame_equal(a.sort_values(key).reset_index(drop=True), b.sort_values(key).reset_index(drop=True))


def test_filterbed_single_chromosome():
    out = run(f"{mawk('filterBed')} {SYNTH}/multi.bed -x -c 21 < {SYNTH}/multi.cpile")
    df = pd.read_csv(StringIO(out), sep="\t")
    assert set(df["Chr"]) == {"chr21"}


@needs_restricted
def test_pon2cnv_legacy_mode_and_first_base_fix(tmp_path):
    snp = tmp_path / "pon.snp"
    cov = run(f"gunzip -c {RESTRICTED}/chr7.pon.gz | {mawk('filterBed')} {DATA}/haev7_padded_chr7.bed -x -L -c chr7 | "
              f"{mawk('PON2CNV')} -x -L -o {snp} -v 0.05 -d 50 -c 0")
    new = pd.read_csv(StringIO(cov), sep="\t")
    old = read(f"{RLEGACY}/pon2cnv_chr7.cov")
    # pooled SNPs: byte-identical
    with open(snp) as a, open(f"{RLEGACY}/pon2cnv_chr7.snp") as b:
        assert a.read() == b.read()
    # same bins; coverage gains the first base of every bin except the very first one
    assert (new[["Chr", "Pos", "ExonPos"]].to_numpy() == old[["Chr", "Pos", "ExonPos"]].to_numpy()).all()
    stream = pd.read_csv(StringIO(run(
        f"gunzip -c {RESTRICTED}/chr7.pon.gz | {mawk('filterBed')} {DATA}/haev7_padded_chr7.bed -x -L -c chr7")), sep="\t")
    depth = stream["Depth"].str.split("=", expand=True)
    d = depth[0].str.split("|", expand=True).astype(int) + depth[1].str.split("|", expand=True).astype(int)
    keep = d.sum(axis=1) > 0
    stream, d = stream.loc[keep].reset_index(drop=True), d.loc[keep].reset_index(drop=True)
    first = (stream["Start"] // 10).drop_duplicates(keep="first").index
    expected = d.loc[first].to_numpy() / 10
    expected[0, :] = 0
    cov_cols = [c for c in new.columns if c.startswith("Cov")]
    delta = new[cov_cols].to_numpy() - old[cov_cols].to_numpy()
    assert np.allclose(delta, expected, atol=1e-9)


def test_cleanpileup_sample_limit_is_numeric(tmp_path):
    # 12 samples; legacy compared -s as a string, so "-s 2" kept all 12
    fields = ["chr1", "100", "A"] + ["3", ".,A", "III"] * 12
    path = tmp_path / "raw.pileup"
    path.write_text("\t".join(fields) + "\n")
    out = run(f"{mawk('cleanpileup')} -s 2 < {path}")
    header = out.splitlines()[0].split("\t")
    assert header == ["Chr", "Start", "Ref", "Read1", "Read2"]


def test_tools_have_bash_shebang():
    from mycnv.tools import mawk_tool

    for name in ["cleanpileup", "filterBed", "pile2CNV", "PON2CNV"]:
        path = mawk_tool(name)
        assert path.read_text().startswith("#!/usr/bin/env bash")
        assert os.access(path, os.X_OK)
