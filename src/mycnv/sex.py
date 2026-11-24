"""Coverage-based sex calls (plan N-Q7, D-Q3; design items M4 and M5).

One detection for the PON columns and for the samples:
- `sample` (and `combine`, `sex`) call the normal and the tumor on their log2ratio vs the (effectively female) PON and
  write `sex.tsv`, whose `ascat_gender` is the gender ASCAT uses. A tumor-only sample (`--tumor-only`) has one row,
  the tumor's, and its call is the ASCAT gender (lower confidence: a chrX loss of a female tumor looks XY).
- `pon-combine` calls every PON normal the same way (`pon_calls`): the log2ratio of its mean-normalised coverage
  (before any doubling) vs the per-bin mean of the other normals, with the XY normals' chrX doubled. The call decides
  which normals get their chrX non-PAR doubled (XY -> XX level, so the PON is effectively female) and is written to
  `CNV/pon.json`.

    log2ratio column: ratio = 2 ** (median(chrX non-PAR) - median(autosomes))   XX ~ 1, XY ~ 0.5
    coverage column:  ratio = median(chrX non-PAR) / median(autosomes)          the same scale, but kit-dependent

XY if ratio < xy_max, XX if ratio > xx_min, else ambiguous; unknown with fewer than min_bins non-PAR chrX bins or no
autosomal bins. NaN and inf values are ignored, and so are zero-coverage bins (a log2ratio is NaN where the coverage
is 0). A bin is on the PAR when Chr is chrX and its Pos (the bin start) lies in a GRCh38 PAR; the other chrX bins are
"chrX non-PAR". Autosomes are chr1..chr22.

A coverage column has no per-bin reference, so its plain ratio is not centred: it carries the kit's chrX/autosome
target mix (HAEv7 PON, chr7 + chr21 + chrX: XX 1.17-1.30, XY 0.61-0.67). `pon_calls` therefore uses the plain ratios
only to find the start of the XY set (`pon_start`: the split into a factor-2 lower and an upper group).
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
# A "split" start of the PON calls needs two groups a factor 2 apart: the log2 distance of the group medians within
# SPLIT_GAP (XY one copy of chrX below XX: 1; HAEv7 PON ~0.97; 47,XXX vs XX 0.585 is outside) and an empty stretch of
# at least SPLIT_STEP between the two groups (HAEv7: XY <= 0.68, XX >= 1.15, i.e. 0.76), so that one outlier and a
# few normals of the main group never pass as the lower group.
SPLIT_GAP = (0.75, 1.25)
SPLIT_STEP = 0.4
CALLS = ("XX", "XY", "ambiguous", "unknown")
SEX_COLUMNS = ["sample", "role", "column", "chrx_ratio", "log2_diff", "n_chrx_bins", "n_auto_bins", "call",
               "ascat_gender"]
_AUTOSOME = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")


def thresholds(cfg: dict | None = None) -> dict:
    """The `sex:` config section (defaults filled in), checked: a mapping, 0 < xy_max <= xx_min, min_bins >= 1."""
    if cfg is not None and not isinstance(cfg, dict):
        raise ValueError(f"the sex: config section must be a mapping (xy_max, xx_min, min_bins), got {cfg!r}")
    t = dict(DEFAULTS)
    t.update(cfg or {})
    try:
        t = {"xy_max": float(t["xy_max"]), "xx_min": float(t["xx_min"]), "min_bins": int(t["min_bins"])}
    except (TypeError, ValueError) as err:
        raise ValueError(f"sex thresholds must be numbers (xy_max, xx_min, min_bins): {err}") from None
    if not 0 < t["xy_max"] <= t["xx_min"]:
        raise ValueError(f"sex thresholds need 0 < xy_max <= xx_min, got xy_max {t['xy_max']}, xx_min {t['xx_min']}")
    if t["min_bins"] < 1:
        raise ValueError(f"sex.min_bins must be >= 1, got {t['min_bins']}")
    return t


def thresholds_record(t: dict, doubling: str | None = None, pon_call: str | None = None,
                      pon_start: str | None = None) -> dict:
    """The thresholds as written to CNV/pon.json (with the PARs and, for the PON, the call, how its XY set started
    and the doubling rule). pon_start "cut" means the PON looked like one sex (see `pon_start`): check the calls."""
    out = dict(t, par_grch38={c: [list(r) for r in regions] for c, regions in PAR.items()})
    if pon_call:
        out["pon_call"] = pon_call
        out["split_gap_log2"] = list(SPLIT_GAP)
        out["split_step_log2"] = SPLIT_STEP
    if pon_start:
        out["pon_start"] = pon_start
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
    """chrx_ratio, log2_diff, n_chrx_bins, n_auto_bins of one column (finite values only; coverage: > 0 only).

    Zero-coverage bins are dropped from a coverage column as the log2ratio path drops them (NaN where Cov is 0): a
    normal with a dropout on one chromosome would otherwise pull its autosomal median down (or to 0).
    """
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    finite = np.isfinite(v) if log2 else np.isfinite(v) & (v > 0)
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


def _values(df: pd.DataFrame, col: str, double: bool, x_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(values with the chrX non-PAR doubled if `double`, usable = finite and > 0) of one column."""
    v = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    if double:
        v = np.where(x_mask, 2 * v, v)
    return v, np.isfinite(v) & (v > 0)


def pon_start(plain: dict, t: dict) -> tuple[set, str]:
    """The XY set the PON calls start from (plain = {column: call_columns row}), and how it was found.

    The kit's chrX target mix scales the plain ratio of every column alike; what sets the XY normals of a mixed PON
    apart is that they sit a factor 2 below the others. So the start is the split of the sorted log2 plain ratios whose
    lower part, raised by 1 (doubled), leaves the PON tightest (least sum of squares around the mean): "split". Only
    splits whose two groups really are a factor 2 apart count: the medians of their log2 plain ratios must differ by
    SPLIT_GAP (0.75-1.25), and no normal may lie within SPLIT_STEP (0.4) of the split. Without that check one normal
    far enough from the others makes "raise everyone else" the tightest split: a 47,XXX normal (1.5x, log2 0.585) in
    an all-XX PON would start, and end, with every XX normal called XY and doubled. A start from the midpoint cut alone fails on a kit with a high chrX mix: an XY ratio >= 0.75
    seeds nothing, and against the undoubled mixed reference a balanced PON's XY normals land in the ambiguous band
    (0.5 / 0.75 ~ 0.67), a fixed point.
    When no split qualifies, the callable normals form one group: the ratios alone cannot tell an all-XX from an all-XY
    PON, and the start is the plain ratio cut at (xy_max + xx_min) / 2, which assumes a chrX/autosome mix near 1 ("cut",
    with a warning). This is right for an all-XX PON on a kit whose plain XX ratio is above the cut and for an all-XY
    PON on a kit whose plain XY ratio is below it; otherwise every normal gets the other sex (documented limit; one
    normal of the other sex in the PON avoids it). Unknown columns never start as XY.
    """
    use = sorted((r["chrx_ratio"], c) for c, r in plain.items()
                 if r["call"] != "unknown" and math.isfinite(r["chrx_ratio"]) and r["chrx_ratio"] > 0)
    n = len(use)
    if n >= 2:
        lg = np.log2([ratio for ratio, _ in use])
        shifted = lg + (np.arange(n)[None, :] < np.arange(n + 1)[:, None])  # row j: the j lowest raised by 1
        ss = ((shifted - shifted.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
        gap = np.array([np.median(lg[j:]) - np.median(lg[:j]) for j in range(1, n)])
        ok = (gap >= SPLIT_GAP[0]) & (gap <= SPLIT_GAP[1]) & (np.diff(lg) >= SPLIT_STEP)
        if ok.any():
            j = int(np.argmin(np.where(ok, ss[1:n], np.inf))) + 1
            if ss[j] < min(ss[0], ss[n]) * (1 - 1e-9):
                return {c for _, c in use[:j]}, "split"
    mid = (t["xy_max"] + t["xx_min"]) / 2
    if n:
        log.warning("PON sex calls: the %d callable normals form one group (plain coverage ratio %.3f-%.3f; no split "
                    "into two groups a factor 2 apart), i.e. one sex: XX or XY is decided by the plain ratio cut at "
                    "%.3f, which depends on the kit's chrX target mix; check the calls against the sheet sex", n,
                    use[0][0], use[-1][0], mid)
    return {c for ratio, c in use if ratio < mid}, "cut"


def pon_calls(df: pd.DataFrame, columns, t: dict, max_rounds: int = 10, start=None) -> list[dict]:
    """Sex calls of PON columns (mean-normalised coverage before any doubling), made like the sample calls.

    Every column is called on its log2ratio vs the per-bin mean of the other PON columns, with the chrX non-PAR of
    the XY columns doubled, i.e. against an effectively female reference (XX ~ 1, XY ~ 0.5, the thresholds of the
    samples). Leaving the column out keeps its call independent of its own doubling. The XY set starts from `start`
    (default: `pon_start`) and is updated until it does not change (at most `max_rounds`). Zero-coverage values enter
    neither the reference nor the column's log2ratio.
    Returns [{column, chrx_ratio, log2_diff, n_chrx_bins, n_auto_bins, call, coverage_ratio, start}]; coverage_ratio
    is the plain ratio (for the log), start how the XY set started ("split", "cut" or "given"; the same in every row).
    """
    columns = list(columns)
    x_mask, auto_mask = region_masks(df)
    plain = {r["column"]: r for r in call_columns(df, columns, t)}
    if start is None:
        xy, how = pon_start(plain, t)
    else:
        xy, how = set(start), "given"
    log.info("PON sex calls: %d columns start as XY (%s)", len(xy), how)
    for rounds in range(1, max_rounds + 1):
        total, count = np.zeros(len(df)), np.zeros(len(df))
        for col in columns:  # per-bin sum and count over all columns (doubled as in `xy`)
            v, ok = _values(df, col, col in xy, x_mask)
            total[ok] += v[ok]
            count[ok] += 1
        calls = []
        for col in columns:
            own, own_ok = _values(df, col, col in xy, x_mask)
            n = count - own_ok
            ref = np.full(len(df), np.nan)
            np.divide(total - np.where(own_ok, own, 0), n, out=ref, where=n > 0)  # mean of the other columns
            v, ok = _values(df, col, False, x_mask)
            ok &= np.isfinite(ref) & (ref > 0)
            lr = np.full(len(df), np.nan)
            lr[ok] = np.log2(v[ok] / ref[ok])
            r = column_ratio(lr, x_mask, auto_mask, log2=True)
            calls.append({"column": col, **r, "call": call_sex(r["chrx_ratio"], r["n_chrx_bins"], r["n_auto_bins"], t),
                          "coverage_ratio": plain[col]["chrx_ratio"], "start": how})
        new = {c["column"] for c in calls if c["call"] == "XY"}
        if new == xy:
            log.info("PON sex calls stable after %d round(s): %d XY", rounds, len(xy))
            break
        xy = new
    else:
        log.warning("PON sex calls still changed after %d rounds: the last calls are used", max_rounds)
    return calls


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
    if cov_df.empty:
        raise ValueError("the coverage table has no bins: no sex call")
    calls = call_columns(cov_df, [col for _, _, col in samples], t, log2=True)
    if not any(r["n_auto_bins"] for r in calls):
        raise ValueError(f"the coverage table has no autosomal bin with a finite "
                         f"{'/'.join(col for _, _, col in samples)}: no sex call")
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
    if "log2ratio1" in cov_df.columns and "log2ratio2" not in cov_df.columns:
        raise ValueError("the coverage table has no log2ratio2 column: one sample column (a tumor-only table? use "
                         "--tumor-only)")
    return sex_table(cov_df, [(normal, "normal", "log2ratio1"), (tumor, "tumor", "log2ratio2")], cfg)


def tumor_sex_table(cov_df: pd.DataFrame, tumor: str, cfg: dict | None = None) -> pd.DataFrame:
    """sex.tsv of a tumor-only sample: one row, tumor = sample 1 (log2ratio1); ascat_gender is the tumor's call.

    Lower confidence than a pair's (the normal's call): a female tumor that lost one X has the chrX level of an XY
    sample and is called XY (a warning says so; the pipeline cross-checks the sheet sex).
    """
    if "log2ratio2" in cov_df.columns:
        raise ValueError("the coverage table has log2ratio1 and log2ratio2 (a T/N pair): tumor-only needs a one-sample "
                         "table (run without --tumor-only)")
    df = sex_table(cov_df, [(tumor, "tumor", "log2ratio1")], cfg)
    if df["call"].iloc[0] == "XY":
        log.warning("Tumor-only sex call XY for %s: a female tumor that lost one X looks the same (lower confidence "
                    "than a normal's call); cross-check with the sheet sex", tumor)
    return df


def json_float(value: float, digits: int = 4):
    """A float for JSON: rounded; NaN/inf -> None (json.dump would write the invalid literal NaN)."""
    return round(float(value), digits) if value is not None and math.isfinite(value) else None
