"""Generic left/right rolling-window engine (main: rollingCNV.py).

For each data column, a trailing window L and a leading window R (L shifted by -(w-1)) are computed, gaps are
interpolated over FullExonPos and the margins are filled from the other side. Then, over the whole table:
Diff = |R - L| / max, value = R * Diff + L * (1 - Diff), optional min-max normalisation, Diff **= diff_exp.

Fixes vs main (CHANGELOG.md):
- windows, interpolation and margin filling run per chromosome (main rolled across chromosome boundaries);
  a chromosome with fewer rows than the window gets partial (trailing L / leading R) windows
- the config key `diff_exp` is read (main read `diffexp` and always used 2; `diffexp` is still accepted)
- `ddof` is passed only to std/var (pandas >= 2 rejects it for mean/sum)
- the debug-mode normalisation of L/R no longer calls a bool
Global steps (Diff scaling, normalisation) stay global, as in main.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .tools import log

_DDOF_AGGS = {"std", "var", "sem"}


def interpolate(df: pd.DataFrame, data_col: str, ref_col: str = "FullExonPos", expand_limit: int = 20) -> pd.DataFrame:
    """Linear interpolation of data_col over the values of ref_col (limit in rows, both directions)."""
    df = df.copy()
    s = pd.Series(df[data_col].to_numpy(), index=pd.Index(df[ref_col].to_numpy()))
    df[data_col] = s.interpolate(method="values", limit=expand_limit, limit_direction="both").to_numpy()
    return df


def normalize_df(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Min-max normalise one column (in a copy)."""
    df = df.copy()
    _min, _max = df[col].min(), df[col].max()
    df[col] = (df[col] - _min) / (_max - _min)
    return df


def llh(data, mean, sigma):
    """Gaussian density of data (Series or array) for the given mean and sigma."""
    s = np.sqrt(2 * np.pi) * sigma
    return np.exp((data - mean) ** 2 / (-2 * (sigma**2))) / s


def _roll(series: pd.Series, aggr, window: int, ddof: int, min_periods: int | None = None) -> pd.Series:
    roll = series.rolling(window, min_periods=min_periods)
    if callable(aggr):
        return roll.apply(aggr, raw=True)
    if aggr in _DDOF_AGGS:
        return roll.agg(aggr, ddof=ddof)
    return roll.agg(aggr)


def _left_right(group: pd.DataFrame, col: str, aggr, window_size: int, ddof: int) -> pd.DataFrame:
    """L and R windows for one chromosome, interpolated and margin-filled within the chromosome."""
    if len(group) >= window_size:
        L = _roll(group[col], aggr, window_size, ddof)
        R = L.shift(-window_size + 1)
    else:
        # chromosome shorter than the window: partial trailing (L) and leading (R) windows
        mp = 2 if callable(aggr) else 1
        L = _roll(group[col], aggr, window_size, ddof, min_periods=mp)
        R = _roll(group[col][::-1], aggr, window_size, ddof, min_periods=mp)[::-1]
    df = pd.DataFrame({"L": L.to_numpy(), "R": R.to_numpy(), "FullExonPos": group["FullExonPos"].to_numpy()})
    for c in ["L", "R"]:
        df = interpolate(df, c, expand_limit=10)
    L, R = df["L"].to_numpy(copy=True), df["R"].to_numpy(copy=True)
    valid_L = np.flatnonzero(~np.isnan(L))
    valid_R = np.flatnonzero(~np.isnan(R))
    if len(valid_L) and len(valid_R):
        first = valid_L[0]
        L[: first + 1] = R[: first + 1]  # main: df.loc[:first_valid_L, "L"] = df["R"] (inclusive)
        last = valid_R[-1]
        R[last + 1 :] = L[last + 1 :]  # main: df.loc[last_valid_R + 1:, "R"] = df["L"]
    else:
        log.warning("Could not fill the rolling margins of %s (%d rows)", col, len(group))
    return pd.DataFrame({"L": L, "R": R}, index=group.index)


def one_col_rolling(df: pd.DataFrame, col: str, aggr="mean", window_size: int = 200, roll_config: dict | None = None,
                    by: str = "Chr") -> pd.DataFrame:
    """Rolling <aggr> of <col>; adds <name> and <name>Diff after <col>.

    aggr is a pandas aggregation name (column name <col>_<aggr>) or a callable taking a numpy window
    (column name = function name). Windows never cross values of the `by` column.
    """
    roll_config = roll_config or {}
    normalize = roll_config.get("normalize", False)
    debug = roll_config.get("debug", False)
    diff_exp = roll_config.get("diff_exp", roll_config.get("diffexp", 2))
    ddof = roll_config.get("ddof", 0)

    cols = list(df.columns)
    col_name = aggr.__name__ if callable(aggr) else f"{col}_{aggr}"
    diff_name = col_name + "Diff"

    parts = [
        _left_right(group, col, aggr, window_size, ddof)
        for _, group in df.groupby(by, sort=False, observed=True)
    ]
    lr = pd.concat(parts).reindex(df.index) if parts else pd.DataFrame({"L": [], "R": []})
    df = df.copy()
    L, R = lr["L"], lr["R"]
    diff = np.abs(R - L)
    diff_max = diff.max()
    diff = diff / diff_max if diff_max > 0 else diff * 0  # main: 0/0 = NaN when L == R everywhere
    df[col_name] = R * diff + L * (1 - diff)
    df[diff_name] = diff
    added = [col_name, diff_name]
    if debug:
        df[f"{col_name}L"], df[f"{col_name}R"] = L, R
        added += [f"{col_name}L", f"{col_name}R"]
    if normalize:
        df = normalize_df(df, col_name)
        if debug:
            for c in [f"{col_name}L", f"{col_name}R"]:
                df = normalize_df(df, c)
    df[diff_name] = df[diff_name] ** diff_exp

    i = cols.index(col) + 1
    return df.loc[:, cols[:i] + added + cols[i:]]


def rolling_data(df: pd.DataFrame, data_params: dict, roll_config: dict | None = None) -> pd.DataFrame:
    """Run one_col_rolling for every column that starts with a key of data_params and every aggregation listed.

    data_params: {column prefix: {aggregation or callable: window size}}. Sum aggregations are min-max
    normalised over the whole table afterwards (as in main).
    """
    roll_config = roll_config or {}
    df = df.copy()
    for data_type, aggs in data_params.items():
        # like main: only the columns present when this data type starts (not the new rolling outputs)
        for data_col in list(df.columns):
            if not data_col.startswith(data_type):
                continue
            for agg, window_size in aggs.items():
                name = agg.__name__ if callable(agg) else agg
                log.info("Rolling %s of %s (window %s)", name, data_col, window_size)
                df = one_col_rolling(df, data_col, aggr=agg, window_size=window_size, roll_config=roll_config)
            # like main: normalise when the last aggregation of this data type is a sum
            if agg != "sum":
                continue
            col_name = f"{data_col}_sum"
            for c in [col_name] + ([f"{col_name}L", f"{col_name}R"] if roll_config.get("debug") else []):
                df = normalize_df(df, c)
    return df
