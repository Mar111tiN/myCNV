"""Step D (PON combine) and step E (sample combine: GC normalisation, FullExonPos, log2ratio vs PON).

Fixes vs main (see CHANGELOG.md):
- PON outlier removal now changes the saved filter stats (stats are recomputed on the outlier-free coverage), and
  pon.cov.full* keep the unmodified coverage. `legacy=True` reproduces the old files exactly.
- FullExonPos offsets follow the natural chromosome order (chr1, chr2, ..., chr22, chrX). main summed them in
  lexicographic order (chr1, chr10, ..., chr19, chr2, ...).
- Chr is a proper ordered categorical (`df["Chr"] = ...`; `df.loc[:, "Chr"] = ...` keeps object dtype on pandas 2+).
- The PON chrX doubling follows the coverage sex call of `mycnv.sex` (XY columns), the same call for the GC and the
  mean table, only on the chrX non-PAR bins (an XY normal has two PAR copies), and is recorded (CNV/pon.json). main
  doubled all chrX bins of columns whose normalised chrX mean was < 75 (`legacy=True`).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import sex
from .ponmeta import check_sample_against_pon, read_pon_meta
from .tools import (
    PRIMARY_CHROMS, chrom_key, chrom_str, log, natural_chroms, read_tsv, set_chrom_categorical, write_tsv,
)


# ##################### PON ###########################################################
def _pon_files(pon_dir: str, sub: str, suffix: str, chroms) -> list[str]:
    files = []
    for chrom in chroms:
        path = os.path.join(pon_dir, sub, f"{chrom}.{suffix}")
        if os.path.isfile(path):
            files.append(path)
        else:
            log.warning("PON file %s not found", path)
    if not files:
        raise FileNotFoundError(f"no {sub}/*.{suffix} files in {pon_dir}")
    return files


def _concat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """pd.concat without the header-only frames: pandas 3 no longer ignores empty frames when it picks the result
    dtypes, so one empty per-chromosome file turned every column into object."""
    full = [d for d in dfs if len(d)]
    return pd.concat(full, ignore_index=True) if full else dfs[0].copy()


def make_pon_snp(pon_dir: str, chroms=PRIMARY_CHROMS) -> pd.DataFrame:
    """Concatenate snp/{chrom}.snp.gz into one table with PONVAF / PONDepth."""
    dfs = [read_tsv(f) for f in _pon_files(pon_dir, "snp", "snp.gz", chroms)]
    snp_df = _concat(dfs).rename(columns={"VAF": "PONVAF", "Depth": "PONDepth"})
    return snp_df


def combine_pon_coverage(pon_dir: str, chroms=PRIMARY_CHROMS) -> pd.DataFrame:
    dfs = [read_tsv(f) for f in _pon_files(pon_dir, "cov", "cov.gz", chroms)]
    cov_df = _concat(dfs)
    if cov_df.empty:
        raise ValueError(f"the PON coverage files in {pon_dir}/cov are empty")
    return set_chrom_categorical(cov_df, chroms)


def normalize_gc_col(cov_df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Scale one coverage column so that every GCratio class has mean 100 on the autosomes (chrX excluded).

    Rows without GCratio (or with a GCratio that only occurs on chrX) are dropped, as in main.
    A class whose autosomal mean is 0 has no defined factor: its bins get NaN (main: 0 * inf = NaN on the autosomes,
    but inf on chrX bins with coverage). This only happens for rare classes, e.g. when few chromosomes are combined.
    """
    norm = cov_df.loc[cov_df["Chr"] != "chrX"].groupby("GCratio", observed=True)[col].mean()
    zero = norm == 0
    if zero.any():
        n_bins = int(cov_df["GCratio"].isin(norm.index[zero]).sum())
        log.warning("GC normalisation of %s: %d GC classes have mean 0 on the autosomes; their %d bins get NaN",
                    col, int(zero.sum()), n_bins)
        norm = norm.mask(zero)
    factor = (100 / norm).rename("factor").reset_index()
    out = cov_df.merge(factor, on="GCratio")
    out[col] = out[col] * out["factor"]
    return out.drop(columns="factor")


def normalize_gc(cov_df: pd.DataFrame) -> pd.DataFrame:
    n = len(cov_df)
    for col in [c for c in cov_df.columns if c.startswith("Cov")]:
        log.info("Normalizing GC ratio for %s", col)
        cov_df = normalize_gc_col(cov_df, col)
    if n and cov_df.empty:
        raise ValueError(f"GC normalisation left none of {n} bins: no bin has a GCratio (missing or wrong GC track)")
    if len(cov_df) < n:
        log.info("GC normalisation drops %d of %d bins without a GCratio of the autosomes", n - len(cov_df), n)
    return cov_df.reset_index(drop=True).sort_values(["Chr", "Pos"])


def normalize(cov_df: pd.DataFrame) -> pd.DataFrame:
    """Scale every coverage column to mean 100 (no GC normalisation)."""
    cov_df = cov_df.copy()
    for col in [c for c in cov_df.columns if c.startswith("Cov")]:
        cov_df[col] = cov_df[col] / cov_df[col].mean() * 100
    return cov_df


def legacy_male_columns(cov_df: pd.DataFrame, x_threshold: float = 75) -> list[str]:
    """main's male rule: Cov columns whose mean over all chrX bins (PAR included) is below `x_threshold`."""
    cov_cols = [c for c in cov_df.columns if c.startswith("Cov")]
    is_x = cov_df["Chr"] == "chrX"
    if not is_x.any():
        return []
    x_mean = cov_df.loc[is_x, cov_cols].mean()
    return [c for c in cov_cols if x_mean[c] < x_threshold]


def amazonize(cov_df: pd.DataFrame, x_threshold: float = 75, columns=None, par: bool = True) -> pd.DataFrame:
    """Double the chrX coverage of male samples, so that they sit at the XX level.

    columns: the columns to double (pon-combine passes the XY calls of `mycnv.sex`); None = main's rule
    (`legacy_male_columns`: normalised chrX mean below `x_threshold`).
    par: also double the PAR bins (main). An XY normal already has two PAR copies, so pon-combine passes False and
    doubles only the chrX non-PAR bins (`mycnv.sex.PAR`).
    """
    cov_df = cov_df.copy()
    is_x = (cov_df["Chr"] == "chrX").to_numpy()
    if not is_x.any():
        return cov_df
    rows = is_x if par else sex.region_masks(cov_df)[0]
    male = legacy_male_columns(cov_df, x_threshold) if columns is None else list(columns)
    if male:
        log.info("Lifting chrX coverage of %d male samples (%s)", len(male), "all chrX bins" if par else
                 f"{int(rows.sum())} chrX non-PAR bins, {int(is_x.sum() - rows.sum())} PAR bins kept")
        cov_df.loc[rows, male] = cov_df.loc[rows, male] * 2
    return cov_df


def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Per-position PONcov_mean / median / std (ddof 1, NaN skipped) over the Cov columns."""
    df = df.drop(columns=[c for c in df.columns if c.startswith("PONcov")])
    cov = df[[c for c in df.columns if c.startswith("Cov")]]
    df = df.copy()
    df["PONcov_mean"] = cov.mean(axis=1)
    df["PONcov_median"] = cov.median(axis=1)
    df["PONcov_std"] = cov.std(axis=1)
    return df


def remove_outliers(df: pd.DataFrame, std_factor: float = 2.5) -> pd.DataFrame:
    """Set coverage values to NaN where |PONcov_mean - value| / PONcov_std > std_factor."""
    df = df.copy()
    for col in [c for c in df.columns if c.startswith("Cov")]:
        df.loc[np.abs(df["PONcov_mean"] - df[col]) / df["PONcov_std"] > std_factor, col] = np.nan
    return df


LEGACY_DOUBLING = "legacy: normalised chrX mean < 75, per table, all chrX bins"
XY_DOUBLING = "XY call, chrX non-PAR bins"
PON_CALL = ("log2ratio vs the per-bin mean of the other PON columns with the chrX non-PAR of the XY columns doubled; "
            "the XY set starts from the split of the sorted log2 coverage ratios whose lower part, doubled, leaves the "
            "PON tightest (one group, i.e. one sex: the coverage ratio cut at (xy_max + xx_min) / 2) and is iterated "
            "until stable")


def pon_sex_calls(cov_mean: pd.DataFrame, sex_cfg: dict | None = None) -> list[dict]:
    """The coverage sex call of every PON column (`mycnv.sex.pon_calls`) on the mean-normalised coverage before
    doubling, against the per-bin mean of the other columns (effectively female: XY chrX doubled)."""
    t = sex.thresholds(sex_cfg)
    cov_cols = [c for c in cov_mean.columns if c.startswith("Cov")]
    calls = sex.pon_calls(cov_mean, cov_cols, t)
    for c in calls:
        plain = c["coverage_ratio"]
        log.info("PON %s: %s (%s; plain coverage ratio %s)", c["column"], c["call"], sex.describe(c),
                 f"{plain:.3f}" if np.isfinite(plain) else "n/a")
    return calls


def make_pon_coverage(pon_dir: str, std_factor: float = 2, chroms=PRIMARY_CHROMS, legacy: bool = False,
                      sex_cfg: dict | None = None) -> dict:
    """Combine cov/{chrom}.cov.gz into the four PON coverage tables (see pon_coverage_and_sex)."""
    return pon_coverage_and_sex(pon_dir, std_factor=std_factor, chroms=chroms, legacy=legacy, sex_cfg=sex_cfg)[0]


def pon_coverage_and_sex(pon_dir: str, std_factor: float = 2, chroms=PRIMARY_CHROMS, legacy: bool = False,
                         sex_cfg: dict | None = None) -> tuple[dict, list[dict]]:
    """Combine cov/{chrom}.cov.gz into the four PON coverage tables, plus the sex call of every PON column.

    Returns ({"full": ..., "fullGC": ..., "filter": ..., "filterGC": ...}, sex calls):
    - full / fullGC: Chr Pos ExonPos Cov1..N PONcov_mean/median/std (mean-normalised / GC-normalised)
    - filter / filterGC: Chr Pos ExonPos + stats after outlier removal
    - sex calls: [{column, chrx_ratio, log2_diff, n_chrx_bins, n_auto_bins, call, coverage_ratio, doubled}]
      (`mycnv.sex.pon_calls`, on the mean-normalised coverage before doubling, vs the effectively female per-bin
      mean of the other columns). The XY columns get their chrX non-PAR bins doubled in both tables (N-Q7; the PAR
      has two copies in XY); XX, ambiguous and unknown columns are not doubled (a warning for each ambiguous or
      unknown one, and one when main's rule would double other columns).
    legacy=True reproduces main: the stats are not recomputed after outlier removal (filter == full stats), the
    full tables carry the outlier-removed (NaN) coverage values, and chrX (all bins) is doubled by main's rule
    (normalised chrX mean < 75, decided separately in the GC and the mean table); the sex calls are then only recorded.
    """
    cov_df = combine_pon_coverage(pon_dir, chroms)
    log.info("Normalizing PON coverage (GC and mean)")
    cov_gc = normalize_gc(cov_df)
    cov_mean = normalize(cov_df)
    calls = pon_sex_calls(cov_mean, sex_cfg)
    xy = [c["column"] for c in calls if c["call"] == "XY"]
    male_gc, male_mean = legacy_male_columns(cov_gc), legacy_male_columns(cov_mean)
    main_rule = set(male_gc) | set(male_mean)
    if set(male_gc) != set(male_mean):
        log.warning("main's chrX rule doubles other columns in the GC table (%s) than in the mean table (%s)",
                    ", ".join(male_gc) or "none", ", ".join(male_mean) or "none")
    if main_rule != set(xy):
        log.warning("main's chrX rule (normalised chrX mean < 75) doubles %s; the coverage sex call has XY for %s%s",
                    ", ".join(sorted(main_rule, key=_cov_key)) or "no column", ", ".join(xy) or "no column",
                    " (--legacy-stats: main's rule is used)" if legacy else " (the XY columns are doubled)")
    if legacy:
        doubled = main_rule
        cov_gc, cov_mean = amazonize(cov_gc, columns=male_gc), amazonize(cov_mean, columns=male_mean)
    else:
        doubled = set(xy)
        for c in calls:
            if c["call"] in ("ambiguous", "unknown"):
                log.warning("PON %s: sex call %s (%s): its chrX is not doubled", c["column"], c["call"],
                            sex.describe(c))
        ambiguous = [c["column"] for c in calls if c["call"] == "ambiguous"]
        if ambiguous:
            log.warning("%d of %d PON columns have an ambiguous sex call and keep their chrX level (%s): check `sex` "
                        "in CNV/pon.json before using this PON", len(ambiguous), len(calls), ", ".join(ambiguous))
        cov_gc = amazonize(cov_gc, columns=xy, par=False)
        cov_mean = amazonize(cov_mean, columns=xy, par=False)
    for c in calls:
        c["doubled"] = c["column"] in doubled
    cov_gc, cov_mean = compute_stats(cov_gc), compute_stats(cov_mean)
    filt_gc, filt_mean = remove_outliers(cov_gc, std_factor), remove_outliers(cov_mean, std_factor)
    if not legacy:
        filt_gc, filt_mean = compute_stats(filt_gc), compute_stats(filt_mean)
    else:
        # main mutated the full tables in place during outlier removal
        cov_gc, cov_mean = filt_gc, filt_mean

    base = ["Chr", "Pos", "ExonPos"]
    cov_cols = [c for c in cov_mean.columns if c.startswith("Cov")]
    stats = ["PONcov_mean", "PONcov_median", "PONcov_std"]
    return {
        "full": chrom_str(cov_mean.loc[:, base + cov_cols + stats]),
        "fullGC": chrom_str(cov_gc.loc[:, base + cov_cols + stats]),
        "filter": chrom_str(filt_mean.loc[:, base + stats]),
        "filterGC": chrom_str(filt_gc.loc[:, base + stats]),
    }, calls


def _cov_key(col: str) -> int:
    return int(col[3:]) if col[3:].isdigit() else 0


def write_pon(pon_dir: str, tables: dict, snp_df: pd.DataFrame) -> list[str]:
    out = []
    for name, df in tables.items():
        path = os.path.join(pon_dir, "CNV", f"pon.cov.{name}.gz")
        log.info("Writing %s", path)
        write_tsv(df, path, compression="gzip")
        out.append(path)
    path = os.path.join(pon_dir, "CNV", "pon.snp.gz")
    write_tsv(chrom_str(snp_df), path, compression="gzip")
    out.append(path)
    return out


# ##################### SAMPLE ########################################################
def concat_tables(paths) -> pd.DataFrame:
    """Concatenate per-chromosome (or whole-sample) tables in natural chromosome order."""
    df = _concat([read_tsv(p) for p in paths])
    return set_chrom_categorical(df).sort_values(["Chr", "Pos"], kind="stable").reset_index(drop=True)


def get_full_exon_pos(*dfs: pd.DataFrame, pon_cov_df: pd.DataFrame, chroms=None):
    """Add FullExonPos = ExonPos + sum of the max PON ExonPos of all preceding chromosomes (natural order)."""
    chrom_max = pon_cov_df.groupby(pon_cov_df["Chr"].astype(str))["ExonPos"].max()
    order = natural_chroms(chrom_max.index) if chroms is None else [c for c in chroms if c in chrom_max.index]
    chrom_add = chrom_max.loc[order].cumsum().shift(fill_value=0).rename("chrAdd")
    out = []
    for df in dfs:
        cols = list(df.columns)
        df = df.copy()
        add = df["Chr"].astype(str).map(chrom_add)
        keep = add.notna()
        if (~keep).any():
            log.warning("Dropping %d rows on chromosomes missing from the PON", int((~keep).sum()))
        df = df.loc[keep].copy()
        df["FullExonPos"] = (df["ExonPos"] + add[keep]).astype(df["ExonPos"].dtype)
        i = cols.index("ExonPos") + 1
        out.append(df.loc[:, cols[:i] + ["FullExonPos"] + cols[i:]].reset_index(drop=True))
    return out


def log2ratio(df: pd.DataFrame, cov_col: str, pon_cov_col: str = "PONcov_mean") -> pd.DataFrame:
    """Insert log2ratio<N> = log2(Cov<N> / PONcov_mean) after Cov<N> (NaN where either value is 0)."""
    cols = list(df.columns)
    df = df.copy()
    log_col = cov_col.replace("Cov", "log2ratio")
    loggable = df[cov_col] * df[pon_cov_col] != 0
    df[log_col] = np.nan
    df.loc[loggable, log_col] = np.log2(df.loc[loggable, cov_col] / df.loc[loggable, pon_cov_col])
    i = cols.index(cov_col) + 1
    return df.loc[:, cols[:i] + [log_col] + cols[i:]]


def _merge_pon_coverage(cov_df: pd.DataFrame, pon_cov_df: pd.DataFrame, legacy_coords: bool,
                        min_pon_match: float) -> pd.DataFrame:
    """Inner-join the sample bins with the PON stats and check that they share one coordinate system.

    Fixed coordinates: bins are keyed by (Chr, Pos) (the bin start, unique per chromosome) and the BED-derived bin
    ExonPos must agree on every shared bin; a difference means another BED or coordinate mode. Legacy coordinates:
    the legacy labels are not unique, so the join stays on (Chr, Pos, ExonPos) as in main.
    In both modes at least `min_pon_match` of the sample bins must find a PON bin.
    """
    n = len(cov_df)
    if legacy_coords:
        merged = cov_df.merge(pon_cov_df, on=["Chr", "Pos", "ExonPos"])
    else:
        dup = pon_cov_df.duplicated(["Chr", "Pos"])
        if dup.any():
            raise ValueError(f"the CNV PON has {int(dup.sum())} duplicate bin labels (Chr, Pos): it was built with "
                             "--legacy-coords (or by the legacy tools); run the sample with --legacy-coords")
        merged = cov_df.merge(pon_cov_df.rename(columns={"ExonPos": "PONExonPos"}), on=["Chr", "Pos"])
        differ = merged["ExonPos"] != merged["PONExonPos"]
        if differ.any():
            first = merged.loc[differ].iloc[0]
            raise ValueError(
                f"{int(differ.sum())} of {len(merged)} bins have another ExonPos in the CNV PON (first "
                f"{first['Chr']}:{first['Pos']}: sample {first['ExonPos']}, PON {first['PONExonPos']}): the sample and "
                "the PON were binned with different BEDs or coordinate modes (--legacy-coords), or the PON was built "
                "by an older mycnv; use the PON's BED and mode or rebuild the CNV PON"
            )
        merged = merged.drop(columns="PONExonPos")
    frac = len(merged) / n if n else 0
    log.info("%d of %d sample bins (%.1f%%) have a PON bin", len(merged), n, 100 * frac)
    if frac < min_pon_match:
        raise ValueError(
            f"only {len(merged)} of {n} sample bins ({100 * frac:.1f}%, minimum {100 * min_pon_match:.0f}%) have a PON "
            "bin: the sample and the PON were binned with different BEDs or coordinate modes (--legacy-coords), or "
            "the PON is for another kit"
        )
    return merged


def combine_sample(cov_df: pd.DataFrame, snp_df: pd.DataFrame, pon_dir: str, gc_normalize: bool = True,
                   legacy_coords: bool = False, min_pon_match: float = 0.9, allow_missing_chroms: bool = False,
                   bed_md5: dict | None = None):
    """Sample coverage/SNPs + PON -> raw tables (main: combine_sample_CNV).

    cov_df: Chr Pos ExonPos GCratio map* Cov1 Cov2 (tn2cnv); snp_df: Chr Pos ExonPos map* Depth1 VAF1 Depth2 VAF2.
    pon_dir: the CNV PON folder (with CNV/). legacy_coords: the sample was binned with --legacy-coords.
    The sample must be complete (every PON chromosome), because GC normalisation and the log2ratio levels are
    relative to the whole sample; `allow_missing_chroms` turns that check off.
    """
    check_sample_against_pon(read_pon_meta(pon_dir), legacy_coords, bed_md5=bed_md5)
    pon_file = os.path.join(pon_dir, "CNV", "pon.cov.filterGC.gz" if gc_normalize else "pon.cov.filter.gz")
    log.info("Loading PON coverage %s", pon_file)
    pon_cov_df = read_tsv(pon_file)
    if pon_cov_df.empty:
        raise ValueError(f"the CNV PON table {pon_file} is empty")
    pon_chroms = set(pon_cov_df["Chr"].astype(str))
    cov_df, snp_df = chrom_str(cov_df), chrom_str(snp_df)
    # main normalised over chr1..chr22, chrX only: bins on other contigs (chrY, alts) must not enter the GC classes
    extra = sorted(set(cov_df["Chr"]) - pon_chroms, key=chrom_key)
    if extra:
        log.info("Dropping %d bins on chromosomes without PON coverage (%s)", int(cov_df["Chr"].isin(extra).sum()),
                 ", ".join(extra))
        cov_df = cov_df.loc[~cov_df["Chr"].isin(extra)]
        snp_df = snp_df.loc[snp_df["Chr"].isin(pon_chroms)]
    missing = sorted(pon_chroms - set(cov_df["Chr"]), key=chrom_key)
    if missing and not allow_missing_chroms:
        raise ValueError(
            f"the sample has no bins on {len(missing)} of {len(pon_chroms)} PON chromosomes ({', '.join(missing[:5])}"
            f"{', ...' if len(missing) > 5 else ''}): GC normalisation and log2ratio need the whole sample, so pass "
            "every per-chromosome tn2cnv output to one `combine` (or --allow-missing-chroms)"
        )
    if not (~cov_df["Chr"].isin(["chrX", "chrY"])).any():
        raise ValueError("the sample has no autosomal bins: GC normalisation uses the autosomes")
    cov_df = set_chrom_categorical(cov_df)
    if gc_normalize:
        log.info("Normalizing sample coverage by GC ratio")
        cov_df = normalize_gc(cov_df)
    else:
        cov_df = normalize(cov_df)
    cov_df, snp_df = chrom_str(cov_df), chrom_str(snp_df)
    cov_df, snp_df = get_full_exon_pos(cov_df, snp_df, pon_cov_df=pon_cov_df)
    cov_df = _merge_pon_coverage(cov_df, pon_cov_df, legacy_coords, min_pon_match)
    del pon_cov_df
    for col in [c for c in cov_df.columns if c.startswith("Cov")]:
        cov_df = log2ratio(cov_df, col)
    pon_snp_file = os.path.join(pon_dir, "CNV", "pon.snp.gz")
    log.info("Merging PON SNPs from %s", pon_snp_file)
    # main merged on all shared columns, i.e. also on the float map* columns; those duplicate the sample's own
    # mappability and break the match when the two files were parsed by different pandas versions
    pon_snp = read_tsv(pon_snp_file)
    pon_snp = pon_snp.drop(columns=[c for c in pon_snp.columns if c.startswith("map")])
    snp_df = snp_df.merge(pon_snp, on=["Chr", "Pos", "ExonPos"], how="left")
    if len(snp_df) and snp_df["PONVAF"].isna().all():
        log.warning("None of the %d sample SNPs is in the PON SNP table %s (PONVAF 0 for all): check that the "
                    "sample and the PON use the same BED and coordinates", len(snp_df), pon_snp_file)
    cov_df = set_chrom_categorical(cov_df).sort_values(["Chr", "Pos"], kind="stable").reset_index(drop=True)
    snp_df = set_chrom_categorical(snp_df).sort_values(["Chr", "Pos"], kind="stable").reset_index(drop=True)
    return chrom_str(cov_df), chrom_str(snp_df)


def filter_cov(cov_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Pre-filter coverage bins for the rolling step (GC range, PON mean and std)."""
    c = config["filter"]["cov"]
    min_gc, max_gc = c["GCrange"]
    return cov_df.loc[
        (cov_df["GCratio"] > min_gc)
        & (cov_df["GCratio"] < max_gc)
        & (cov_df["PONcov_mean"] >= c["minPONcov"])
        & (cov_df["PONcov_std"] < c["maxPONstd"])
    ]
