"""Step F: centre-Gaussian likelihood of the log2ratios and rolling coverage windows (main: rollingCov.py)."""

from __future__ import annotations

import pandas as pd

from .rolling import interpolate, llh, rolling_data
from .tools import log, set_chrom_categorical, chrom_str


def compute_covLLH_col(df: pd.DataFrame, log_col: str, params: dict) -> pd.DataFrame:
    """covLLH<N>: density of log2ratio<N> under the Gaussian fitted to the centre band (whole sample)."""
    cols = list(df.columns)
    lo, hi = params["center_range"]
    center = df.loc[(df[log_col] > lo) & (df[log_col] < hi), log_col]
    sigma = center.std() * params["sigma_factor"]
    mean = center.mean()
    log.info("covLLH of %s: centre gaussian mean %.3f sigma %.3f", log_col, mean, sigma)
    llh_col = log_col.replace("log2ratio", "covLLH")
    df = df.copy()
    df[llh_col] = llh(df[log_col], mean=mean, sigma=sigma)
    i = cols.index(log_col) + 1
    return df.loc[:, cols[:i] + [llh_col] + cols[i:]]


def compute_covLLH(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    for col in list(df.columns):
        if col.startswith("log2ratio"):
            df = compute_covLLH_col(df, col, params=config["rolling"]["cov"]["LLH"])
    return df


def rolling_coverage(cov_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Rolling windows over the filtered coverage bins, per chromosome in (Chr, Pos) order."""
    cc = config["rolling"]["cov"]
    data_params = cc["data"]
    df = set_chrom_categorical(cov_df).sort_values(["Chr", "Pos"], kind="stable")
    if "covLLH" in data_params:
        df = compute_covLLH(df, config)
    df = rolling_data(df, data_params=data_params, roll_config=cc)
    return chrom_str(df)


def interpolate_fullexonpos(merge_df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing FullExonPos per chromosome by interpolation over Pos (main: interpolate_fullexonpon)."""
    parts = [
        interpolate(chrom_df, "FullExonPos", ref_col="Pos", expand_limit=1_000_000)
        for _, chrom_df in merge_df.groupby("Chr", sort=False, observed=True)
    ]
    df = pd.concat(parts).sort_values("FullExonPos", kind="stable")
    df["FullExonPos"] = df["FullExonPos"].fillna(0).astype(int)
    return df
