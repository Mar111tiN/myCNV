"""Step D (PON combine) and step E (sample combine: GC normalisation, FullExonPos, log2ratio vs PON).

Fixes vs main (see CHANGELOG.md):
- PON outlier removal now changes the saved filter stats (stats are recomputed on the outlier-free coverage), and
  pon.cov.full* keep the unmodified coverage. `legacy=True` reproduces the old files exactly.
- FullExonPos offsets follow the natural chromosome order (chr1, chr2, ..., chr22, chrX). main summed them in
  lexicographic order (chr1, chr10, ..., chr19, chr2, ...).
- Chr is a proper ordered categorical (`df["Chr"] = ...`; `df.loc[:, "Chr"] = ...` keeps object dtype on pandas 2+).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

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


def make_pon_snp(pon_dir: str, chroms=PRIMARY_CHROMS) -> pd.DataFrame:
    """Concatenate snp/{chrom}.snp.gz into one table with PONVAF / PONDepth."""
    dfs = [read_tsv(f) for f in _pon_files(pon_dir, "snp", "snp.gz", chroms)]
    snp_df = pd.concat(dfs, ignore_index=True).rename(columns={"VAF": "PONVAF", "Depth": "PONDepth"})
    return snp_df


def combine_pon_coverage(pon_dir: str, chroms=PRIMARY_CHROMS) -> pd.DataFrame:
    dfs = [read_tsv(f) for f in _pon_files(pon_dir, "cov", "cov.gz", chroms)]
    return set_chrom_categorical(pd.concat(dfs, ignore_index=True), chroms)


def normalize_gc_col(cov_df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Scale one coverage column so that every GCratio class has mean 100 on the autosomes (chrX excluded).

    Rows without GCratio (or with a GCratio that only occurs on chrX) are dropped, as in main.
    """
    norm = cov_df.loc[cov_df["Chr"] != "chrX"].groupby("GCratio", observed=True)[col].mean()
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


def amazonize(cov_df: pd.DataFrame, x_threshold: float = 75) -> pd.DataFrame:
    """Double the chrX coverage of samples whose chrX mean is below `x_threshold` (male samples -> XX level)."""
    cov_df = cov_df.copy()
    cov_cols = [c for c in cov_df.columns if c.startswith("Cov")]
    is_x = cov_df["Chr"] == "chrX"
    if not is_x.any():
        return cov_df
    x_mean = cov_df.loc[is_x, cov_cols].mean()
    male = [c for c in cov_cols if x_mean[c] < x_threshold]
    if male:
        log.info("Lifting chrX coverage of %d male samples", len(male))
        cov_df.loc[is_x, male] = cov_df.loc[is_x, male] * 2
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


def make_pon_coverage(pon_dir: str, std_factor: float = 2, chroms=PRIMARY_CHROMS, legacy: bool = False) -> dict:
    """Combine cov/{chrom}.cov.gz into the four PON coverage tables.

    Returns {"full": ..., "fullGC": ..., "filter": ..., "filterGC": ...}:
    - full / fullGC: Chr Pos ExonPos Cov1..N PONcov_mean/median/std (mean-normalised / GC-normalised)
    - filter / filterGC: Chr Pos ExonPos + stats after outlier removal
    legacy=True reproduces main: the stats are not recomputed after outlier removal (filter == full stats) and
    the full tables carry the outlier-removed (NaN) coverage values.
    """
    cov_df = combine_pon_coverage(pon_dir, chroms)
    log.info("Normalizing PON coverage (GC and mean)")
    cov_gc = normalize_gc(cov_df)
    cov_mean = normalize(cov_df)
    cov_gc, cov_mean = amazonize(cov_gc), amazonize(cov_mean)
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
    }


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
    df = pd.concat([read_tsv(p) for p in paths], ignore_index=True)
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
