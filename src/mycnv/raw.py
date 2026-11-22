"""Step B/C: pileup (sample) or PONmatrix (PON) -> 10 bp coverage bins + het-SNP / pooled-SNP tables.

Both wrap the mawk chain `filterBed | pile2CNV` (sample) or `filterBed | PON2CNV` (PON) and then merge the GC ratio
and the mappability tracks.
"""

from __future__ import annotations

import os
import tempfile

import pandas as pd

from .tools import cat_cmd, cmd2df, log, mawk
from .tracks import add_gc, add_genmap


def _filter_cmd(bed: str, chrom: str | None, legacy_coords: bool) -> str:
    cmd = f"{mawk('filterBed')} {bed} -x"
    if chrom:
        cmd += f" -c {chrom}"
    if legacy_coords:
        cmd += " -L"
    return cmd


def tn2cnv(
    pileup: str,
    bed: str,
    gc: str | None = None,
    genmap: str | None = None,
    chrom: str | None = None,
    normal_vaf=(0.2, 0.8),
    min_depth: int = 15,
    min_cov: float = 0,
    bin_size: int = 10,
    clean: bool = True,
    legacy_coords: bool = False,
    tmp_dir: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage bins and het SNPs of a paired pileup (normal = sample 1).

    pileup: `samtools mpileup ... NORMAL TUMOR | cleanpileup -d` output (clean=True, the somVarWES pileup files),
            or raw `samtools mpileup` text (clean=False; cleanpileup -d is added). Plain or gzipped.
    chrom:  restrict to one chromosome (optional split, plan D16); default is the whole sample.
    legacy_coords: legacy bin labels and ExonPos (only to match PON coverage built by the legacy tools).
    """
    read_cmd = cat_cmd(pileup)
    if not clean:
        read_cmd += f" | {mawk('cleanpileup')} -d"
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        snp_file = os.path.join(tmp, "sample.snp")
        cnv_cmd = (
            f"{mawk('pile2CNV')} -x -o {snp_file} -v {normal_vaf[0]} -V {normal_vaf[1]} "
            f"-d {min_depth} -c {min_cov} -b {bin_size}"
        )
        if legacy_coords:
            cnv_cmd += " -L"
        cmd = f"{read_cmd} | {_filter_cmd(bed, chrom, legacy_coords)} | {cnv_cmd}"
        log.info("Coverage and het SNPs from %s%s", pileup, f" ({chrom})" if chrom else "")
        cov_df = cmd2df(cmd)
        snp_df = pd.read_csv(snp_file, sep="\t")
    cov_df = add_gc(cov_df, gc)
    cov_df, snp_df = add_genmap(cov_df, snp_df, track=genmap)
    return cov_df, snp_df


def pon2cnv(
    matrix: str,
    bed: str,
    gc: str | None = None,
    genmap: str | None = None,
    chrom: str | None = None,
    min_vaf: float = 0.05,
    min_depth: int = 50,
    min_cov: float = 0,
    bin_size: int = 10,
    legacy_coords: bool = False,
    tmp_dir: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample coverage bins and pooled SNP VAFs of one PONmatrix file (`matrix/{chrom}.pon.gz`)."""
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        snp_file = os.path.join(tmp, "pon.snp")
        cnv_cmd = f"{mawk('PON2CNV')} -x -o {snp_file} -v {min_vaf} -d {min_depth} -c {min_cov} -b {bin_size}"
        if legacy_coords:
            cnv_cmd += " -L"
        cmd = f"{cat_cmd(matrix)} | {_filter_cmd(bed, chrom, legacy_coords)} | {cnv_cmd}"
        log.info("PON coverage and SNPs from %s", matrix)
        cov_df = cmd2df(cmd)
        snp_df = pd.read_csv(snp_file, sep="\t")
    cov_df = add_gc(cov_df, gc)
    cov_df, snp_df = add_genmap(cov_df, snp_df, track=genmap)
    return cov_df, snp_df
