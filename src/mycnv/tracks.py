"""GC-ratio and mappability (genmap) tracks: load per chromosome and merge onto coverage / SNP tables.

A track is given either as a directory of per-chromosome files (the legacy "split" layout, file name pattern with
`{chrom}`) or as a single file: a full-genome file with a `Chr` column, or one chromosome's file (no `Chr` column; its
name must name the chromosome, e.g. `chr21.gc100-10.gz`).

Missing data fail loudly (main crashed; the first port turned this into empty outputs with exit 0): a track path
that does not exist, no data for a primary chromosome (chr1-22, chrX) that has rows, or a track that matches none of
a primary chromosome's rows raise an error. Contigs outside chr1-22/X (chrY, alts) only get a warning; their bins are
dropped later (the CNV PON has no such chromosomes).
"""

from __future__ import annotations

import os
import re
from functools import lru_cache

import pandas as pd

from .tools import log

GC_PATTERN = "{chrom}.gc100-10.gz"
GENMAP_PATTERN = "hg38_genmap.{chrom}.txt.gz"
GENMAP_MODES = ["30_0", "50_0", "75_1", "100_2"]
PRIMARY_RE = re.compile(r"^chr([0-9]+|X)$")
MIN_ROWS_ALL_NAN = 10  # a primary-chromosome table with at least this many rows and no track value is an error


def is_primary(chrom: str) -> bool:
    return bool(PRIMARY_RE.match(str(chrom)))


def check_track(track: str | None, kind: str) -> str:
    if not track:
        raise ValueError(f"no {kind} track given")
    if not os.path.exists(track):
        raise FileNotFoundError(f"{kind} track {track} does not exist")
    return track


def _names_chrom(path: str, chrom: str) -> bool:
    return bool(re.search(rf"(?<![A-Za-z0-9]){re.escape(chrom)}(?![A-Za-z0-9])", os.path.basename(path)))


def _no_data(kind: str, chrom: str, track: str) -> None:
    msg = f"no {kind} data for {chrom} in {track}"
    if is_primary(chrom):
        raise FileNotFoundError(f"{msg} (wrong track, file name pattern or chromosome naming?)")
    log.warning("%s: its rows get no %s value", msg, kind)


@lru_cache(maxsize=4)
def _read_full(path: str, usecols: tuple | None) -> pd.DataFrame:
    log.info("Loading track %s", path)
    return pd.read_csv(path, sep="\t", usecols=list(usecols) if usecols else None)


def _usecols(path: str, usecols: tuple | None) -> tuple | None:
    if usecols is None:
        return None
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    return tuple(c for c in usecols if c in header)


def _load(track: str, chrom: str, pattern: str, kind: str, usecols: tuple | None) -> pd.DataFrame | None:
    """Rows of one chromosome from a track directory or file; None if there are none."""
    if os.path.isdir(track):
        path = os.path.join(track, pattern.format(chrom=chrom))
        if not os.path.isfile(path):
            return None
        cols = _usecols(path, usecols)
        df = pd.read_csv(path, sep="\t", usecols=list(cols) if cols else None)
    else:
        df = _read_full(track, _usecols(track, usecols))
        if "Chr" in df.columns:
            df = df.loc[df["Chr"].astype(str) == chrom]
        elif not _names_chrom(track, chrom):
            raise ValueError(f"{kind} file {track} has no Chr column and its name does not name {chrom}: give a "
                             "directory of per-chromosome files, a file with a Chr column, or the chromosome's file")
    return None if df.empty else df


def load_gc(track: str, chrom: str, pattern: str = GC_PATTERN) -> pd.DataFrame | None:
    """GC table for one chromosome with columns [Chr,] Pos, GCratio (from `Start`, `GC/AT`); None if no data."""
    df = _load(track, chrom, pattern, "GC", None)
    if df is None:
        return None
    df = df.rename(columns={"GC/AT": "GCratio", "Start": "Pos"})
    return df.loc[:, [c for c in ["Chr", "Pos", "GCratio"] if c in df.columns]]


def load_genmap(track: str, chrom: str, modes=GENMAP_MODES, pattern: str = GENMAP_PATTERN) -> pd.DataFrame | None:
    """Mappability for one chromosome: Chr, Pos, map<mode>...; gaps forward- then back-filled (as legacy)."""
    cols = ["Chr", "Pos"] + list(modes)
    df = _load(track, chrom, pattern, "mappability", tuple(cols))
    if df is None:
        return None
    if "Chr" not in df.columns:
        df = df.assign(Chr=chrom)
    df = df.loc[:, cols].ffill().bfill()
    return df.rename(columns={m: f"map{m}" for m in modes})


def _check_merged(df: pd.DataFrame, chrom: str, col: str, kind: str, track: str, what: str) -> None:
    n_nan = int(df[col].isna().sum())
    if not n_nan:
        return
    if n_nan == len(df) and len(df) >= MIN_ROWS_ALL_NAN and is_primary(chrom):
        raise ValueError(f"the {kind} track {track} has no value at any of the {len(df)} {what} on {chrom}: wrong "
                         "track (other chromosome, genome build or kit split?)")
    log.warning("%d of %d %s on %s have no %s value (%s)", n_nan, len(df), what, chrom, kind,
                "dropped by the GC normalisation" if kind == "GC" else "gaps in the track; the SNP filter drops them")


def add_gc(cov_df: pd.DataFrame, track: str | None, pattern: str = GC_PATTERN) -> pd.DataFrame:
    """Left-merge GCratio onto a coverage table, chromosome by chromosome."""
    track = check_track(track, "GC")
    parts = []
    for chrom, df in cov_df.groupby("Chr", sort=False, observed=True):
        chrom = str(chrom)
        gc_df = load_gc(track, chrom, pattern)
        if gc_df is None:
            _no_data("GC", chrom, track)
            df = df.assign(GCratio=float("nan"))
        else:
            df = df.merge(gc_df, on=["Chr", "Pos"] if "Chr" in gc_df.columns else "Pos", how="left")
            _check_merged(df, chrom, "GCratio", "GC", track, "bins")
        parts.append(df)
    return pd.concat(parts, ignore_index=True) if parts else cov_df.assign(GCratio=pd.Series(dtype=float))


def add_genmap(*dfs: pd.DataFrame, track: str | None, modes=GENMAP_MODES, pattern: str = GENMAP_PATTERN,
               required: bool = True, names=("bins", "SNPs")):
    """Left-merge mappability onto each table and move the base and map columns to the front (legacy order).

    required=False (PON coverage: the map columns are not used downstream) leaves the map columns empty when no
    track is given.
    """
    map_cols = [f"map{m}" for m in modes]
    if not track and not required:
        log.info("No mappability track: map columns stay empty (not used by the CNV PON)")
    else:
        track = check_track(track, "mappability")
    out = []
    for k, df in enumerate(dfs):
        what = names[k] if k < len(names) else "rows"
        if not track:
            merged = df.assign(**{c: float("nan") for c in map_cols})
        else:
            parts = []
            for chrom, part in df.groupby("Chr", sort=False, observed=True):
                chrom = str(chrom)
                gm = load_genmap(track, chrom, modes, pattern)
                if gm is None:
                    _no_data("mappability", chrom, track)
                    part = part.assign(**{c: float("nan") for c in map_cols})
                else:
                    part = part.merge(gm, on=["Chr", "Pos"], how="left")
                    _check_merged(part, chrom, map_cols[0], "mappability", track, what)
                parts.append(part)
            merged = pd.concat(parts, ignore_index=True) if parts else df.assign(**{c: pd.Series(dtype=float) for c in map_cols})
        cols = list(merged.columns)
        base = [c for c in cols if c in ["Chr", "Start", "Pos", "ExonPos", "GCratio"] or c.startswith("map")]
        data = [c for c in cols if c not in base]
        out.append(merged.loc[:, base + data])
    return out
