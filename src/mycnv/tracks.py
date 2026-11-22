"""GC-ratio and mappability (genmap) tracks: load per chromosome and merge onto coverage / SNP tables.

A track is given either as a directory of per-chromosome files (the legacy "split" layout, file name pattern with
`{chrom}`) or as a single full-genome file with a `Chr` column.
"""

from __future__ import annotations

import os
from functools import lru_cache

import pandas as pd

from .tools import log

GC_PATTERN = "{chrom}.gc100-10.gz"
GENMAP_PATTERN = "hg38_genmap.{chrom}.txt.gz"
GENMAP_MODES = ["30_0", "50_0", "75_1", "100_2"]


def _track_file(track: str, chrom: str, pattern: str) -> str | None:
    if os.path.isdir(track):
        path = os.path.join(track, pattern.format(chrom=chrom))
        return path if os.path.isfile(path) else None
    return track if os.path.isfile(track) else None


@lru_cache(maxsize=4)
def _read_full(path: str, usecols: tuple | None) -> pd.DataFrame:
    log.info("Loading full track %s", path)
    return pd.read_csv(path, sep="\t", usecols=list(usecols) if usecols else None)


def load_gc(track: str, chrom: str, pattern: str = GC_PATTERN) -> pd.DataFrame | None:
    """GC table for one chromosome with columns [Chr,] Pos, GCratio (from `Start`, `GC/AT`)."""
    path = _track_file(track, chrom, pattern)
    if path is None:
        log.warning("No GC data for %s in %s", chrom, track)
        return None
    if os.path.isdir(track):
        df = pd.read_csv(path, sep="\t")
    else:
        df = _read_full(path, None)
        df = df.loc[df["Chr"] == chrom] if "Chr" in df.columns else df
    return df.rename(columns={"GC/AT": "GCratio", "Start": "Pos"})


def load_genmap(track: str, chrom: str, modes=GENMAP_MODES, pattern: str = GENMAP_PATTERN) -> pd.DataFrame | None:
    """Mappability for one chromosome: Chr, Pos, map<mode>...; gaps forward- then back-filled (as legacy)."""
    path = _track_file(track, chrom, pattern)
    if path is None:
        log.warning("No mappability data for %s in %s", chrom, track)
        return None
    cols = ["Chr", "Pos"] + list(modes)
    if os.path.isdir(track):
        df = pd.read_csv(path, sep="\t", usecols=cols)
    else:
        df = _read_full(path, tuple(cols))
        df = df.loc[df["Chr"] == chrom]
    df = df.loc[:, cols].ffill().bfill()
    return df.rename(columns={m: f"map{m}" for m in modes})


def add_gc(cov_df: pd.DataFrame, track: str | None, pattern: str = GC_PATTERN) -> pd.DataFrame:
    """Left-merge GCratio onto a coverage table, chromosome by chromosome."""
    if not track:
        log.warning("No GC track given: GCratio stays empty")
        out = cov_df.copy()
        out["GCratio"] = float("nan")
        return out
    parts = []
    for chrom, df in cov_df.groupby("Chr", sort=False, observed=True):
        gc_df = load_gc(track, str(chrom), pattern)
        if gc_df is None:
            df = df.assign(GCratio=float("nan"))
        elif "Chr" in gc_df.columns:
            df = df.merge(gc_df, on=["Chr", "Pos"], how="left")
        else:
            df = df.merge(gc_df, on="Pos", how="left")
        parts.append(df)
    return pd.concat(parts, ignore_index=True) if parts else cov_df.assign(GCratio=pd.Series(dtype=float))


def add_genmap(*dfs: pd.DataFrame, track: str | None, modes=GENMAP_MODES, pattern: str = GENMAP_PATTERN):
    """Left-merge mappability onto each table and move the base and map columns to the front (legacy order)."""
    out = []
    map_cols = [f"map{m}" for m in modes]
    for df in dfs:
        if not track:
            log.warning("No mappability track given: map columns stay empty")
            merged = df.assign(**{c: float("nan") for c in map_cols})
        else:
            parts = []
            for chrom, part in df.groupby("Chr", sort=False, observed=True):
                gm = load_genmap(track, str(chrom), modes, pattern)
                if gm is None:
                    part = part.assign(**{c: float("nan") for c in map_cols})
                else:
                    part = part.merge(gm, on=["Chr", "Pos"], how="left")
                parts.append(part)
            merged = pd.concat(parts, ignore_index=True) if parts else df.assign(**{c: pd.Series(dtype=float) for c in map_cols})
        cols = list(merged.columns)
        base = [c for c in cols if c in ["Chr", "Start", "Pos", "ExonPos", "GCratio"] or c.startswith("map")]
        data = [c for c in cols if c not in base]
        out.append(merged.loc[:, base + data])
    return out
