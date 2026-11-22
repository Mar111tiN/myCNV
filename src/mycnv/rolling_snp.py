"""Step G: merge rolled coverage onto het SNPs, centre the VAFs, SNP density / fallSNP, SNP filter, remerge.

(main: rollingSNP.py). Per-chromosome work (interpolation, rolling) never crosses chromosomes; the global steps
(VAF centring, normalisations, fallSNP scaling) stay global. Rows are processed in (Chr, FullExonPos) order.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .rolling import interpolate, rolling_data
from .tools import chrom_str, log, set_chrom_categorical

LOG2_PAT = re.compile(r"log2ratio[0-9]+(_mean)?$")


def _sort(df: pd.DataFrame) -> pd.DataFrame:
    return set_chrom_categorical(df).sort_values(["Chr", "FullExonPos"], kind="stable").reset_index(drop=True)


def merge_cov_to_snp(snp_df: pd.DataFrame, cov_df: pd.DataFrame) -> pd.DataFrame:
    """Outer-merge the rolled coverage onto the SNPs, interpolate the coverage columns per chromosome over
    FullExonPos (limit 200 rows) and keep the SNP rows (main: mergeCov2SNP)."""
    snp_df = chrom_str(snp_df)
    cov_df = chrom_str(cov_df)
    base_cols = list(snp_df.columns)[:4]  # Chr Pos ExonPos FullExonPos
    cov_cols = [c for c in cov_df.columns if c.startswith("PON") or LOG2_PAT.match(c)] + ["GCratio"]
    cov_df = cov_df.loc[:, base_cols + cov_cols]
    snp_df = snp_df.copy()
    for col in ["PONVAF", "PONDepth"]:
        snp_df[col] = snp_df[col].fillna(0)
    df = _sort(snp_df.merge(cov_df, how="outer"))
    parts = []
    for _, chrom_df in df.groupby("Chr", sort=False, observed=True):
        for col in cov_cols:
            chrom_df = interpolate(chrom_df, col, expand_limit=200)
        parts.append(chrom_df)
    df = pd.concat(parts) if parts else df
    return df.loc[df["PONVAF"].notna()].reset_index(drop=True)


def center_vaf(snp_df: pd.DataFrame) -> pd.DataFrame:
    """Shift every VAF<N> so that its mean is 0.5; add absVAF<N> and offVAF<N> = 2 * (VAF - 0.5)."""
    df = snp_df.copy()
    for col in [c for c in df.columns if c.startswith("VAF")]:
        cols = list(df.columns)
        i = cols.index(col) + 1
        df[col] = df[col] + 0.5 - df[col].mean()
        df[f"off{col}"] = (df[col] - 0.5) * 2
        df[f"abs{col}"] = np.abs(df[f"off{col}"])
        df = df.loc[:, cols[:i] + [f"abs{col}", f"off{col}"] + cols[i:]]
    return df


def make_get_density(window_size: int = 20):
    """Rolling callback: SNPs per exonic base in the window (window_size / span of FullExonPos)."""

    def SNPdensity(data):
        return window_size / (data.max() - data.min())

    return SNPdensity


def get_fall_snp(snp_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """SNPdensity (rolling, normalised) and fallSNP = SNPdensity * PONVAF / max."""
    c = config["rolling"]["snp"]
    window = c["data"]["SNPdensity"]
    params = {"FullExonPos": {make_get_density(window): window}}
    df = rolling_data(snp_df, data_params=params, roll_config=c)
    df["fallSNP"] = df["SNPdensity"] * df["PONVAF"]
    top = df["fallSNP"].max()
    if top > 0:
        df["fallSNP"] = df["fallSNP"] / top
    else:
        # main divided by 0 (NaN everywhere), so the fallSNP filter dropped every SNP; PONVAF 0 = no PON artefact
        log.warning("No SNP has a PON VAF: fallSNP set to 0 (are the sample and the PON on the same BED and "
                    "coordinates?)")
        df["fallSNP"] = 0.0
    return df


def filter_snp(snp_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Keep SNPs with PONVAF < maxPONVAF, mappability >= thresholds, low fallSNP and SNP density."""
    c = config["filter"]["snp"]
    df = snp_df.copy()
    for col in ["PONVAF", "PONDepth"]:
        df[col] = df[col].fillna(0)
    keep = (df["PONVAF"] < c["maxPONVAF"]) & (df["fallSNP"] < c["maxFallSNP"]) & (df["SNPdensity"] < c["maxSNPdensity"])
    map_keys = [k for k in c if k.startswith("map")]
    no_map = df[map_keys].isna().any(axis=1) if map_keys else pd.Series(False, index=df.index)
    if no_map.any():
        log.warning("%d of %d SNPs have no mappability value (gaps in the genmap track): the SNP filter drops them",
                    int(no_map.sum()), len(df))
    for key in map_keys:
        keep &= df[key] >= c[key]
    log.info("SNP filter keeps %d of %d SNPs", int(keep.sum()), len(df))
    return df.loc[keep].reset_index(drop=True)


def rolling_snp(snp_df: pd.DataFrame, roll_cov_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    log.info("Merging SNP and coverage data")
    df = merge_cov_to_snp(snp_df, roll_cov_df)
    log.info("Centering VAFs")
    df = center_vaf(df)
    log.info("SNP density and fallSNP")
    df = get_fall_snp(df, config)
    df = filter_snp(df, config)
    return chrom_str(df)


def remerge_cnv(snp_df: pd.DataFrame, cov_df: pd.DataFrame) -> pd.DataFrame:
    """SNP rows and coverage rows in one track (coverage rows have empty VAFs) (main: remergeCNV)."""
    base_cols = ["Chr", "Pos", "ExonPos", "FullExonPos"]
    snp_cols = [c for c in snp_df.columns if c.startswith("VAF")]
    cov_cols = [c for c in cov_df.columns if LOG2_PAT.match(c)]
    cnv_df = pd.concat(
        [chrom_str(snp_df).loc[:, base_cols + snp_cols + cov_cols], chrom_str(cov_df).loc[:, base_cols + cov_cols]]
    )
    return chrom_str(_sort(cnv_df))
