"""Coverage-based sex calls (plan N-Q7, D-Q3; design items M4 and M5).

One detection for the PON columns and for the samples:
- `pon-combine` calls every PON normal on its mean-normalised coverage (before any doubling). The call decides which
  normals get their chrX doubled (XY -> XX level, so the PON is effectively female) and is written to `CNV/pon.json`.
- `sample` (and `combine --out-sex`, `sex`) call the normal and the tumor on their log2ratio vs that PON and write
  `sex.tsv`, whose `ascat_gender` is the gender ASCAT uses.

    coverage column:  ratio = median(chrX non-PAR) / median(autosomes)          XX ~ 1, XY ~ 0.5
    log2ratio column: ratio = 2 ** (median(chrX non-PAR) - median(autosomes))   the same scale, the same thresholds

XY if ratio < xy_max, XX if ratio > xx_min, else ambiguous; unknown with fewer than min_bins non-PAR chrX bins or no
autosomal bins. NaN and inf values are ignored. A bin is on the PAR when Chr is chrX and its Pos (the bin start) lies
in a GRCh38 PAR; the other chrX bins are "chrX non-PAR". Autosomes are chr1..chr22.
"""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from .tools import log

# GRCh38 pseudo-autosomal regions on chrX, 1-based inclusive (XY: two copies, like an autosome)
PAR = {"chrX": ((10_001, 2_781_479), (155_701_383, 156_030_895))}
DEFAULTS = {"xy_max": 0.65, "xx_min": 0.85, "min_bins": 100}
CALLS = ("XX", "XY", "ambiguous", "unknown")
SEX_COLUMNS = ["sample", "role", "column", "chrx_ratio", "log2_diff", "n_chrx_bins", "n_auto_bins", "call",
               "ascat_gender"]
_AUTOSOME = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")


def thresholds(cfg: dict | None = None) -> dict:
    """The `sex:` config section (defaults filled in), checked: 0 < xy_max <= xx_min, min_bins >= 1."""
    t = dict(DEFAULTS)
    t.update(cfg or {})
    t = {"xy_max": float(t["xy_max"]), "xx_min": float(t["xx_min"]), "min_bins": int(t["min_bins"])}
    if not 0 < t["xy_max"] <= t["xx_min"]:
        raise ValueError(f"sex thresholds need 0 < xy_max <= xx_min, got xy_max {t['xy_max']}, xx_min {t['xx_min']}")
    if t["min_bins"] < 1:
        raise ValueError(f"sex.min_bins must be >= 1, got {t['min_bins']}")
    return t


def thresholds_record(t: dict, doubling: str | None = None) -> dict:
    """The thresholds as written to CNV/pon.json (with the PARs and, for the PON, the doubling rule)."""
    out = dict(t, par_grch38={c: [list(r) for r in regions] for c, regions in PAR.items()})
    if doubling:
        out["doubling"] = doubling
    return out


def region_masks(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """(chrX non-PAR, autosome) boolean masks over the rows of a table with Chr and Pos."""
    chrom = df["Chr"].astype(str).to_numpy()
    pos = pd.to_numeric(df["Pos"], errors="coerce").to_numpy(float)
    is_x = chrom == "chrX"
    par = np.zeros(len(df), bool)
    for start, end in PAR["chrX"]:
        par |= is_x & (pos >= start) & (pos <= end)
    autosomes = [c for c in pd.unique(chrom) if _AUTOSOME.match(c)]
    return is_x & ~par, np.isin(chrom, autosomes)


def column_ratio(values, x_mask: np.ndarray, auto_mask: np.ndarray, log2: bool = False) -> dict:
    """chrx_ratio, log2_diff, n_chrx_bins, n_auto_bins of one column (finite values only)."""
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    finite = np.isfinite(v)
    x, a = v[x_mask & finite], v[auto_mask & finite]
    ratio = diff = math.nan
    if len(x) and len(a):
        if log2:
            diff = float(np.median(x) - np.median(a))
            ratio = float(2 ** diff)
        else:
            mx, ma = float(np.median(x)), float(np.median(a))
            if ma > 0:
                ratio = mx / ma
                diff = math.log2(ratio) if ratio > 0 else -math.inf
    return {"chrx_ratio": ratio, "log2_diff": diff, "n_chrx_bins": int(len(x)), "n_auto_bins": int(len(a))}


def call_sex(ratio: float, n_chrx_bins: int, n_auto_bins: int, t: dict) -> str:
    if n_chrx_bins < t["min_bins"] or n_auto_bins == 0 or not math.isfinite(ratio):
        return "unknown"
    if ratio < t["xy_max"]:
        return "XY"
    if ratio > t["xx_min"]:
        return "XX"
    return "ambiguous"


def call_columns(df: pd.DataFrame, columns, t: dict, log2: bool = False) -> list[dict]:
    """[{column, chrx_ratio, log2_diff, n_chrx_bins, n_auto_bins, call}] for the given columns of df."""
    x_mask, auto_mask = region_masks(df)
    out = []
    for col in columns:
        r = column_ratio(df[col], x_mask, auto_mask, log2=log2)
        out.append({"column": col, **r, "call": call_sex(r["chrx_ratio"], r["n_chrx_bins"], r["n_auto_bins"], t)})
    return out


def describe(r: dict) -> str:
    ratio = "n/a" if not math.isfinite(r["chrx_ratio"]) else f"{r['chrx_ratio']:.3f}"
    return f"chrX/autosome ratio {ratio}, {r['n_chrx_bins']} non-PAR chrX bins, {r['n_auto_bins']} autosomal bins"


def ascat_gender(rows: list[dict]) -> str:
    """The gender for ASCAT: the normal's call (pairs) or the tumor's (no normal row); ambiguous/unknown -> XX."""
    if not rows:
        raise ValueError("no sample rows for the sex call")
    ref = next((r for r in rows if r["role"] == "normal"), rows[0])
    if ref["call"] in ("XX", "XY"):
        return ref["call"]
    log.warning("Sex call of the %s %s is %s (%s): ASCAT gender XX", ref["role"], ref["sample"], ref["call"],
                describe(ref))
    return "XX"


def sex_table(cov_df: pd.DataFrame, samples, cfg: dict | None = None) -> pd.DataFrame:
    """sex.tsv rows for samples = [(label, role, log2ratio column), ...] of a raw coverage table (log2ratio vs PON)."""
    t = thresholds(cfg)
    missing = [col for _, _, col in samples if col not in cov_df.columns]
    if missing:
        raise ValueError(f"the coverage table has no {', '.join(missing)} column: the sex call needs the raw coverage "
                         "table of `combine`/`sample` (raw/<s>.cov.gz, log2ratio vs the PON)")
    calls = call_columns(cov_df, [col for _, _, col in samples], t, log2=True)
    rows = []
    for (label, role, _), r in zip(samples, calls):
        rows.append({"sample": label, "role": role, **r})
        log.info("Sex call of the %s %s (%s): %s (%s)", role, label, r["column"], r["call"], describe(r))
    gender = ascat_gender(rows)
    for row in rows:
        row["ascat_gender"] = gender
    df = pd.DataFrame(rows).loc[:, SEX_COLUMNS]
    df[["chrx_ratio", "log2_diff"]] = df[["chrx_ratio", "log2_diff"]].astype(float).round(4)
    return df


def pair_sex_table(cov_df: pd.DataFrame, tumor: str, normal: str, cfg: dict | None = None) -> pd.DataFrame:
    """sex.tsv of a T/N pair: normal = sample 1 (log2ratio1), tumor = sample 2 (log2ratio2)."""
    return sex_table(cov_df, [(normal, "normal", "log2ratio1"), (tumor, "tumor", "log2ratio2")], cfg)


def json_float(value: float, digits: int = 4):
    """A float for JSON: rounded; NaN/inf -> None (json.dump would write the invalid literal NaN)."""
    return round(float(value), digits) if value is not None and math.isfinite(value) else None
