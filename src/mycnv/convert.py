"""Step H: ASCAT input files (BAF and logR of normal and tumor, plus a mock germline) (main: convert.write_ASCAT).

Files (tab-separated, pandas index as unnamed first column, which ASCAT reads as row names):
  {outdir}/{sample}_{normal}_baf.tsv   {sample}_{tumor}_baf.tsv   {sample}_germline_baf.tsv
  {outdir}/{sample}_{normal}_logr.tsv  {sample}_{tumor}_logr.tsv  {sample}_germline_logr.tsv
Columns: chrs (without "chr"), pos, <sample>. BAF = centred VAF<N>; logR = log2ratio<N>_mean (rolling-smoothed).
"""

from __future__ import annotations

import os
import re

import pandas as pd

from .tools import log

LOG_MEAN_PAT = re.compile(r"log2ratio[0-9]+_mean$")


def split_sample_name(sample: str) -> tuple[str, str]:
    """Legacy naming <patient>_<tumor>-<normal> (e.g. NHL1_A-B) -> (tumor, normal)."""
    try:
        tumor, normal = sample.rsplit("_", 1)[1].split("-")
    except (ValueError, IndexError) as err:
        raise ValueError(f"cannot derive tumor/normal from {sample!r}; pass them explicitly") from err
    return tumor, normal


def write_ascat(cnv_df: pd.DataFrame, sample: str, outdir: str, tumor: str | None = None, normal: str | None = None,
                mock_germline: bool = True) -> list[str]:
    """Write the ASCAT BAF/logR tables; normal = sample 1 (VAF1, log2ratio1_mean), tumor = sample 2."""
    if tumor is None or normal is None:
        tumor, normal = split_sample_name(sample)
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, sample)
    df = cnv_df.copy()
    df["Chr"] = df["Chr"].astype(str).str.replace("chr", "")
    vaf_cols = [c for c in df.columns if c.startswith("VAF")]
    log_cols = [c for c in df.columns if LOG_MEAN_PAT.match(c)]
    if len(vaf_cols) < 2 or len(log_cols) < 2:
        raise ValueError(f"need VAF1/VAF2 and log2ratio1_mean/log2ratio2_mean, got {vaf_cols} {log_cols}")

    def _write(col: str, path: str, values=None) -> str:
        out = df.loc[:, ["Chr", "Pos", col]].rename(columns={"Chr": "chrs", "Pos": "pos", col: sample})
        if values is not None:
            out[sample] = values
        out.to_csv(path, sep="\t")
        return path

    files = [
        _write(vaf_cols[0], f"{base}_{normal}_baf.tsv"),
        _write(vaf_cols[1], f"{base}_{tumor}_baf.tsv"),
    ]
    if mock_germline:
        files.append(_write(vaf_cols[0], f"{base}_germline_baf.tsv", values=0.5))
    files += [
        _write(log_cols[0], f"{base}_{normal}_logr.tsv"),
        _write(log_cols[1], f"{base}_{tumor}_logr.tsv"),
    ]
    if mock_germline:
        files.append(_write(log_cols[0], f"{base}_germline_logr.tsv", values=0))
    log.info("ASCAT input written to %s_*", base)
    return files
