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
    """Legacy naming <patient>_<tumor>-<normal> (e.g. NHL1_A-B) -> (tumor, normal).

    Only the legacy form is parsed: the last "_" field must be <tumor>-<normal>. Pair IDs like sarek's
    <tumor>_vs_<normal> parsed into wrong labels without an error (P01-T_vs_P01-N -> P01, N), so they are rejected:
    pass the names explicitly.
    """
    err = f"cannot derive tumor/normal from {sample!r} (legacy form <patient>_<tumor>-<normal>); pass " \
          "--tumor-name and --normal-name"
    if not sample or "_vs_" in sample:
        raise ValueError(err)
    try:
        tumor, normal = sample.rsplit("_", 1)[1].split("-")
    except (ValueError, IndexError) as e:
        raise ValueError(err) from e
    if not tumor or not normal:
        raise ValueError(err)
    return tumor, normal


def resolve_names(sample: str | None, tumor: str | None, normal: str | None) -> tuple[str, str]:
    """Both names given -> them; none -> parsed from the legacy sample ID; one alone is an error."""
    if (tumor is None) != (normal is None):
        raise ValueError("give both --tumor-name and --normal-name, or neither (then they come from the sample ID)")
    if tumor is None:
        tumor, normal = split_sample_name(sample)
        log.info("Tumor/normal names from the sample ID %s: %s / %s", sample, tumor, normal)
    return tumor, normal


def write_ascat(cnv_df: pd.DataFrame, sample: str, outdir: str, tumor: str | None = None, normal: str | None = None,
                mock_germline: bool = True, allow_empty: bool = False) -> list[str]:
    """Write the ASCAT BAF/logR tables; normal = sample 1 (VAF1, log2ratio1_mean), tumor = sample 2.

    Raises ValueError when there is no SNP row or when a logR column has no value (ASCAT cannot segment that; main
    wrote such tables and ASCAT failed later), unless allow_empty is set.
    """
    tumor, normal = resolve_names(sample, tumor, normal)
    if tumor == normal or (mock_germline and "germline" in (tumor, normal)):
        raise ValueError(f"tumor and normal names {tumor!r}/{normal!r} would give clashing ASCAT file names")
    base = os.path.join(outdir, sample)
    df = cnv_df.copy()
    df["Chr"] = df["Chr"].astype(str).str.replace("chr", "")
    vaf_cols = [c for c in df.columns if c.startswith("VAF")]
    log_cols = [c for c in df.columns if LOG_MEAN_PAT.match(c)]
    if len(vaf_cols) < 2 or len(log_cols) < 2:
        raise ValueError(f"need VAF1/VAF2 and log2ratio1_mean/log2ratio2_mean, got {vaf_cols} {log_cols}")
    problems = []
    if df.empty:
        problems.append("no SNP passed the SNP filter")
    for col in log_cols[:2]:
        n_nan = int(df[col].isna().sum())
        if len(df) and n_nan == len(df):
            problems.append(f"{col} (logR) has no value")
        elif n_nan:
            log.warning("%d of %d SNP rows have no %s (logR)", n_nan, len(df), col)
    if problems:
        msg = f"ASCAT input of {sample}: " + "; ".join(problems)
        if not allow_empty:
            raise ValueError(msg + " (use --allow-empty to write the tables anyway)")
        log.warning("%s: writing the tables anyway (--allow-empty)", msg)
    os.makedirs(outdir, exist_ok=True)

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
