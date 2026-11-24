"""Candidate het-SNP sites for tumor-only samples (`--snp-sites`, plan D27 / N-Q3, design item M1).

A tumor-only sample has no normal to say which positions are germline het SNPs, so the candidates come from a
common-SNP sites file (the kit's `<kit>.snp_sites.vcf.gz`: biallelic SNVs from gnomAD / dbSNP common sites, cut to the
kit's padded BED). Two formats are read, told apart by their content (not by the file name), plain or gzipped/bgzipped:

- VCF: only CHROM, POS, REF and ALT are read. Biallelic SNVs are kept (REF and ALT one base A/C/G/T); multi-allelic
  records (ALT with a comma), indels and symbolic alleles are dropped. A position listed twice (a split multi-allelic
  site) counts once.
- TSV with the two columns Chr and Pos (1-based, as in the pileup), with or without a header line `Chr Pos`.

Chromosome names are matched as written (no "chr" prefix is added or removed). The file is read in chunks and only
the positions that are candidates in the sample are kept, so a large (even genome-wide) sites file is fine.
"""

from __future__ import annotations

import csv
import gzip

import numpy as np
import pandas as pd

from .tools import chrom_key, log

CHUNK = 500_000  # rows per read: bounds the memory of a large sites file
_SNV = r"[ACGTacgt]"


def _is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _sniff(path: str, gz: bool) -> tuple[str, int, list[str] | None]:
    """(format "vcf" or "tsv", number of lines before the first data line, fields of the first data line)."""
    is_vcf, skip, first = False, 0, None
    with (gzip.open(path, "rt") if gz else open(path)) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                is_vcf |= line.startswith(("##fileformat=VCF", "#CHROM"))
                skip += 1
                continue
            first = line.rstrip("\r\n").split("\t")
            break
    return ("vcf" if is_vcf else "tsv"), skip, first


def _chunks(path: str):
    """(format, iterator over (CHROM, POS, REF, ALT) column chunks as strings; REF/ALT None for a TSV)."""
    gz = _is_gzip(path)
    fmt, skip, first = _sniff(path, gz)
    kw = dict(sep="\t", header=None, quoting=csv.QUOTE_NONE, na_filter=False, chunksize=CHUNK,
              compression="gzip" if gz else None)
    if first is None:
        return fmt, iter(())
    if fmt == "vcf":
        if len(first) < 5:
            raise ValueError(f"--snp-sites {path}: a VCF data line needs CHROM POS ID REF ALT, got {len(first)} fields")
        reader = pd.read_csv(path, skiprows=skip, usecols=[0, 1, 3, 4], dtype={0: str, 1: str, 3: str, 4: str}, **kw)
        return fmt, ((c[0], c[1], c[3], c[4]) for c in reader)
    if len(first) != 2:
        raise ValueError(f"--snp-sites {path}: not a VCF (no ##fileformat=VCF / #CHROM line), and a sites TSV has the "
                         f"2 columns Chr Pos (1-based), got {len(first)} (a BED is not accepted: it is 0-based)")
    if not first[1].strip().isdigit():  # a header line such as "Chr Pos"
        skip += 1
    reader = pd.read_csv(path, skiprows=skip, usecols=[0, 1], dtype={0: str, 1: str}, **kw)
    return fmt, ((c[0], c[1], None, None) for c in reader)


def _positions(pos: pd.Series, path: str) -> np.ndarray:
    num = pd.to_numeric(pos, errors="coerce")
    bad = num.isna() | (num < 1) | (num != np.floor(num))
    if bad.any():
        raise ValueError(f"--snp-sites {path}: position {pos[bad].iloc[0]!r} is not a 1-based integer")
    return num.to_numpy("int64")


def sites_at(path: str, candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """The sites of `path` that are among `candidates` (Chr, Pos), and read statistics.

    Returns (DataFrame Chr, Pos of the matching sites, unique), {"format", "records", "sites" (kept by the format
    rules), "chroms" (set of the site chromosomes)}.
    """
    cand = {str(c): np.unique(g["Pos"].to_numpy("int64")) for c, g in candidates.groupby("Chr", observed=True)}
    fmt, chunks = _chunks(path)
    stats = {"format": fmt, "records": 0, "sites": 0, "chroms": set()}
    found = []
    for chrom, pos, ref, alt in chunks:
        stats["records"] += len(chrom)
        if fmt == "vcf":
            snv = (ref.str.fullmatch(_SNV) & alt.str.fullmatch(_SNV)).to_numpy(bool)
            chrom, pos = chrom[snv], pos[snv]
        stats["sites"] += len(chrom)
        positions = _positions(pos, path)
        chroms = chrom.to_numpy(str)
        for c in pd.unique(chroms):
            stats["chroms"].add(c)
            if c not in cand:
                continue
            p = positions[chroms == c]
            hit = p[np.isin(p, cand[c])]
            if len(hit):
                found.append(pd.DataFrame({"Chr": c, "Pos": hit}))
    out = pd.concat(found, ignore_index=True).drop_duplicates() if found else pd.DataFrame(
        {"Chr": pd.Series(dtype=str), "Pos": pd.Series(dtype="int64")})
    return out.reset_index(drop=True), stats


def _shown(names) -> str:
    names = sorted(names, key=chrom_key)
    return ", ".join(names[:5]) + (", ..." if len(names) > 5 else "")


def restrict_to_sites(snp_df: pd.DataFrame, path: str, chroms) -> pd.DataFrame:
    """Keep the rows of a het-SNP candidate table (Chr, Pos, ...) whose position is in the sites file `path`.

    Fails (ValueError) when the file has no usable site, or when none of its chromosomes is one of `chroms` (the
    BED's: naming, e.g. "7" vs "chr7", or a sites file of another genome build).
    """
    matched, st = sites_at(path, snp_df[["Chr", "Pos"]])
    what = "biallelic SNV records" if st["format"] == "vcf" else "sites"
    log.info("SNP sites %s (%s): %d records, %d %s, %d at a candidate position", path, st["format"].upper(),
             st["records"], st["sites"], what, len(matched))
    if not st["sites"]:
        raise ValueError(f"--snp-sites {path} has no {what} ({st['records']} records): wrong file?")
    chroms = {str(c) for c in chroms}
    if not st["chroms"] & chroms:
        raise ValueError(f"--snp-sites {path} has no site on the chromosomes of the BED ({_shown(chroms)}; sites on "
                         f"{_shown(st['chroms'])}): chromosome naming or genome build of the sites file?")
    keys = pd.MultiIndex.from_frame(matched.astype({"Chr": str}))
    keep = pd.MultiIndex.from_frame(snp_df[["Chr", "Pos"]].astype({"Chr": str, "Pos": "int64"})).isin(keys)
    return snp_df.loc[keep].reset_index(drop=True)
