"""Provenance of a CNV PON: coordinate mode, bin size and BED regions, so samples cannot be mixed with a PON binned
another way (review 2026-11-22: a default sample against a --legacy-coords PON merged 0 bins and exited 0).

- `pon-coverage` writes `cov/{chrom}.cov.json` next to every `cov/{chrom}.cov.gz` (one file per chromosome, so
  per-chromosome tasks can be gathered into one folder without name clashes).
- `pon-combine` checks that they agree and writes `CNV/pon.json`.
- `combine`/`sample` read `CNV/pon.json` (if present) and fail on a coordinate-mode or BED mismatch. PONs built by the
  legacy tools have no metadata; then only the data checks in `combine.combine_sample` apply.
"""

from __future__ import annotations

import json
import os

from . import __version__
from .tools import log

COORDS = {False: "fixed", True: "legacy"}


def cov_meta_path(out_dir: str, chrom: str) -> str:
    return os.path.join(out_dir, "cov", f"{chrom}.cov.json")


def write_cov_meta(out_dir: str, chrom: str, legacy_coords: bool, bin_size: int, bed: str, bed_md5: str) -> str:
    path = cov_meta_path(out_dir, chrom)
    meta = {"mycnv": __version__, "chrom": chrom, "coords": COORDS[bool(legacy_coords)], "bin_size": bin_size,
            "bed": os.path.basename(bed), "bed_regions_md5": bed_md5}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=1)
    return path


def combine_cov_meta(pon_dir: str, chroms) -> dict | None:
    """Merge cov/{chrom}.cov.json of the given chromosomes; None if none exists (legacy PON)."""
    metas = {}
    for chrom in chroms:
        path = cov_meta_path(pon_dir, chrom)
        if os.path.isfile(path):
            with open(path) as fh:
                metas[chrom] = json.load(fh)
    if not metas:
        log.warning("No cov/*.cov.json in %s (PON coverage built by the legacy tools?): the coordinate mode of the "
                    "CNV PON is not recorded", pon_dir)
        return None
    for key in ["coords", "bin_size"]:
        values = {m[key] for m in metas.values()}
        if len(values) > 1:
            raise ValueError(f"PON coverage files in {pon_dir}/cov disagree on {key}: {sorted(map(str, values))}")
    first = next(iter(metas.values()))
    return {"coords": first["coords"], "bin_size": first["bin_size"], "bed": sorted({m["bed"] for m in metas.values()}),
            "bed_regions_md5": {c: m["bed_regions_md5"] for c, m in metas.items()},
            "coverage_mycnv": sorted({m["mycnv"] for m in metas.values()})}


def write_pon_meta(out_dir: str, meta: dict | None, chroms, std_factor: float, legacy_stats: bool) -> str:
    path = os.path.join(out_dir, "CNV", "pon.json")
    full = {"mycnv": __version__, "chroms": list(chroms), "std_factor": std_factor, "legacy_stats": legacy_stats}
    full.update(meta or {"coords": None})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(full, fh, indent=1)
    return path


def read_pon_meta(pon_dir: str) -> dict | None:
    path = os.path.join(pon_dir, "CNV", "pon.json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def check_sample_against_pon(meta: dict | None, legacy_coords: bool, bin_size: int = 10,
                             bed_md5: dict | None = None) -> None:
    """Fail when the sample was binned differently from the PON (coordinate mode, bin size, BED regions)."""
    if not meta or meta.get("coords") is None:
        return
    mode = COORDS[bool(legacy_coords)]
    if meta["coords"] != mode:
        hint = "add --legacy-coords" if meta["coords"] == "legacy" else "drop --legacy-coords"
        raise ValueError(f"the CNV PON uses {meta['coords']} coordinates, the sample {mode} ones: {hint} (or rebuild "
                         "the CNV PON from matrix/ in the sample's mode)")
    if meta.get("bin_size") not in (None, bin_size):
        raise ValueError(f"the CNV PON uses {meta['bin_size']} bp bins, the sample {bin_size} bp")
    pon_md5 = meta.get("bed_regions_md5") or {}
    if bed_md5:
        differ = sorted(c for c in set(pon_md5) & set(bed_md5) if pon_md5[c] != bed_md5[c])
        if differ:
            raise ValueError(f"the sample BED differs from the BED the CNV PON was built with ({', '.join(meta.get('bed', []))}) "
                             f"on {', '.join(differ)}: use the same BED for the PON and the samples")
