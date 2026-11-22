"""Capture BED regions as `filterBed.mawk` reads them, and the BED-derived bin ExonPos of the fixed coordinates.

filterBed (fixed coordinates) gives the stream position `pos` of region i the ExonPos
    BEDSUM[i] - BEDEND[i] + pos,      BEDSUM[i] = sum over regions k <= i of the chromosome of (end_k - start_k + 1),
i.e. a region [start, end) covers the stream positions start..end (both inclusive, kept from v1.0.1). A position is
assigned to the first region (in BED order) whose end is >= pos, and is kept if pos >= that region's start.

pile2CNV/PON2CNV label a bin by its start (Pos = 10 * floor(pos / 10)) but give it the ExonPos of the first base that
the stream holds in the bin. That depends on the data (PON2CNV skips bases without coverage in every normal, a pileup
restricted by `samtools mpileup -l` lacks the base at the 0-based BED start, a tumor can have reads where the PON has
none), so the same bin could get different ExonPos values in a sample and in the PON. `bin_exon_pos` replaces it with
the ExonPos of the bin's first BED base, which depends on the BED alone.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

Regions = dict  # {chrom: (starts, ends)} as int64 arrays, in BED order


def normalize_chrom(chrom: str) -> str:
    """filterBed `-c` convention: "7" -> "chr7"."""
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def read_bed(path: str, chroms=None) -> Regions:
    """{chrom: (starts, ends)} of a BED file, read with the rules of filterBed.mawk.

    Lines starting with browser/track/# and empty lines are skipped; fields are tab-separated; only the first three
    columns are used. filterBed needs every chromosome as one block sorted by start (the blocks may come in any
    order); otherwise it silently drops or misassigns positions, so this raises ValueError instead.
    """
    keep = None if chroms is None else {normalize_chrom(c) for c in chroms}
    starts: dict[str, list[int]] = {}
    ends: dict[str, list[int]] = {}
    seen: set[str] = set()
    last_chrom, last_start = None, None
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith(("browser", "track", "#")) or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                raise ValueError(f"{path}:{n}: not a tab-separated BED line: {line.rstrip()!r}")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom != last_chrom:
                if chrom in seen:
                    raise ValueError(f"{path}:{n}: {chrom} comes back after {last_chrom}: each chromosome must be one "
                                     "block (sort the BED with `sort -k1,1 -k2,2n`)")
                seen.add(chrom)
            elif start < last_start:
                raise ValueError(f"{path}:{n}: {chrom}:{start} is listed after {chrom}:{last_start}: sort the BED "
                                 "by start within each chromosome (`sort -k1,1 -k2,2n`)")
            last_chrom, last_start = chrom, start
            if keep is not None and chrom not in keep:
                continue
            starts.setdefault(chrom, []).append(start)
            ends.setdefault(chrom, []).append(end)
    return {c: (np.asarray(starts[c], dtype=np.int64), np.asarray(ends[c], dtype=np.int64)) for c in starts}


def region_digest(regions: Regions, chrom: str) -> str:
    """md5 of one chromosome's regions (start, end pairs in BED order): equal digests = same bins and ExonPos."""
    s, e = regions[chrom]
    text = "".join(f"{a}\t{b}\n" for a, b in zip(s.tolist(), e.tolist()))
    return hashlib.md5(text.encode()).hexdigest()


def bin_exon_pos(df: pd.DataFrame, regions: Regions, bin_size: int = 10) -> np.ndarray:
    """ExonPos of the first BED base of every bin (fixed coordinates: `Pos` = bin start).

    df needs Chr and Pos. Raises ValueError for a bin that holds no BED base (the BED is not the one the bins were
    made with).
    """
    out = np.zeros(len(df), dtype=np.int64)
    pos_all = df["Pos"].to_numpy(dtype=np.int64)
    for chrom, idx in df.groupby(df["Chr"].astype(str), sort=False).indices.items():
        if chrom not in regions:
            raise ValueError(f"bins on {chrom}, but the BED has no {chrom} region")
        s, e = regions[chrom]
        cum = np.cumsum(e - s + 1)
        first_end = np.maximum.accumulate(e)  # region of a position = first region (BED order) with end >= pos
        b0 = pos_all[idx]
        i = np.searchsorted(first_end, b0, side="left")
        bad = i >= len(s)
        i = np.minimum(i, len(s) - 1)
        p0 = np.maximum(b0, s[i])
        bad |= p0 >= b0 + bin_size
        if bad.any():
            raise ValueError(f"{int(bad.sum())} bins on {chrom} (first at {int(b0[np.argmax(bad)])}) hold no base of "
                             "the BED: the bins were made with another BED")
        out[idx] = cum[i] - e[i] + p0
    return out
