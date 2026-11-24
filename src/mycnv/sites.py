"""Candidate het-SNP sites for tumor-only samples (`--snp-sites`, plan D27 / N-Q3, design item M1).

A tumor-only sample has no normal to say which positions are germline het SNPs, so the candidates come from a
common-SNP sites file (the kit's `<kit>.snp_sites.vcf.gz`: biallelic SNVs from gnomAD / dbSNP common sites, cut to the
kit's padded BED). Two formats are read, told apart by their content (not by the file name), plain or gzipped/bgzipped:

- VCF: a file with a `##fileformat=VCF` line or a `#CHROM POS ID REF ALT ...` column header (5 or more fields). Only
  CHROM, POS, REF and ALT are read. Biallelic SNVs are kept (REF and ALT one base A/C/G/T); indels and symbolic
  alleles are dropped. A multi-allelic site is dropped however it is written: as one record (ALT with a comma) or
  split into biallelic records at one position (`bcftools norm -m-`, gnomAD v4), i.e. a position with records of
  more than one allele is dropped. An identical record listed twice counts once.
- TSV with the two columns Chr and Pos (1-based, as in the pileup). `#` lines are comments (so `#CHROM POS` is a
  header); a first line whose Pos field has no digit (`Chr Pos`) is the header; any other line is data, and a Pos
  that is not a 1-based integer is an error.

Chromosome names are matched as written (no "chr" prefix is added or removed). The file is read in chunks and only
the positions that are candidates in the sample are kept, so a large (even genome-wide) sites file is fine. It must be
a regular file (it is read more than once).
"""

from __future__ import annotations

import csv
import gzip
import os

import numpy as np
import pandas as pd

from .tools import chrom_key, log

CHUNK = 500_000  # rows per read: bounds the memory of a large sites file
_SNV = r"[ACGTacgt]"


def _is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _is_vcf_header(line: str) -> bool:
    """`##fileformat=VCF...` or the VCF column header `#CHROM POS ID REF ALT ...` (a 2-column `#CHROM POS`: a TSV)."""
    return line.startswith("##fileformat=VCF") or (
        line.startswith("#CHROM") and len(line.rstrip("\r\n").split("\t")) >= 5)


def _sniff(path: str, gz: bool) -> tuple[str, int, list[str] | None]:
    """(format "vcf" or "tsv", number of lines before the first data line, fields of the first data line or None).

    Comment (`#`) and blank lines and a TSV header (the first other line, 2 fields, a Pos field without a digit) are
    no data lines: a header-only file has none (None).
    """
    is_vcf, skip, first, header_allowed = False, 0, None, True
    with (gzip.open(path, "rt") if gz else open(path)) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                is_vcf |= _is_vcf_header(line)
                skip += 1
                continue
            fields = line.rstrip("\r\n").split("\t")
            if header_allowed and not is_vcf and len(fields) == 2 and not any(ch.isdigit() for ch in fields[1]):
                header_allowed = False  # only the first line can be the header: a second "Chr Pos" line is an error
                skip += 1
                continue
            first = fields
            break
    return ("vcf" if is_vcf else "tsv"), skip, first


def _layout(path: str) -> tuple[bool, str, int, list[str] | None]:
    """(gzipped, format, lines before the first data line, fields of the first data line or None); ValueError when
    the first data line is no VCF record (a VCF) or has not the 2 columns Chr Pos (otherwise)."""
    gz = _is_gzip(path)
    fmt, skip, first = _sniff(path, gz)
    if first is not None and fmt == "vcf" and len(first) < 5:
        raise ValueError(f"--snp-sites {path}: a VCF data line needs CHROM POS ID REF ALT, got {len(first)} fields")
    if first is not None and fmt == "tsv" and len(first) != 2:
        raise ValueError(f"--snp-sites {path}: not a VCF (no ##fileformat=VCF line or #CHROM POS ID REF ALT header), "
                         f"and a sites TSV has the 2 columns Chr Pos (1-based), got {len(first)} (a BED is not "
                         "accepted: it is 0-based)")
    return gz, fmt, skip, first


def check_sites(path: str) -> str:
    """Check a sites file before the pileup pass: it exists, is a regular file, has the layout of a VCF or a Chr Pos
    TSV and at least one data line (FileNotFoundError / ValueError otherwise). Returns the format."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"--snp-sites {path} does not exist")
    if not os.path.isfile(path):
        raise ValueError(f"--snp-sites {path} is not a regular file (it is read more than once: no pipe, FIFO or "
                         "process substitution)")
    _, fmt, _, first = _layout(path)
    if first is None:
        raise ValueError(f"--snp-sites {path} has no records (empty or header only): wrong file?")
    return fmt


def _chunks(path: str):
    """(format, iterator over (CHROM, POS, REF, ALT) column chunks as strings; REF/ALT None for a TSV)."""
    gz, fmt, skip, first = _layout(path)
    if first is None:
        return fmt, iter(())
    kw = dict(sep="\t", header=None, skiprows=skip, quoting=csv.QUOTE_NONE, na_filter=False, chunksize=CHUNK,
              compression="gzip" if gz else None)
    if fmt == "vcf":
        reader = pd.read_csv(path, usecols=[0, 1, 3, 4], dtype={0: str, 1: str, 3: str, 4: str}, **kw)
        return fmt, ((c[0], c[1], c[3], c[4]) for c in reader)
    reader = pd.read_csv(path, usecols=[0, 1], dtype={0: str, 1: str}, **kw)
    return fmt, ((c[0], c[1], None, None) for c in reader)


def _positions(pos: pd.Series, path: str, check: np.ndarray) -> np.ndarray:
    """Positions as int64. A position that is not a 1-based integer is an error in the rows of `check` (the kept
    records) and -1 (never a candidate) in the others (e.g. a VCF telomere record at POS 0)."""
    num = pd.to_numeric(pos, errors="coerce").to_numpy("float64")
    bad = np.isnan(num) | (num < 1) | (num != np.floor(num))
    if (bad & check).any():
        raise ValueError(f"--snp-sites {path}: position {pos.iloc[np.flatnonzero(bad & check)[0]]!r} is not a 1-based "
                         "integer")
    return np.where(bad, -1, num).astype("int64")


def sites_at(path: str, candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """The sites of `path` that are among `candidates` (Chr, Pos), and read statistics.

    Returns (DataFrame Chr, Pos of the matching sites, unique), {"format", "records", "sites" (records kept by the
    format rules: biallelic SNVs of a VCF, every row of a TSV), "chroms" (set of their chromosomes), "multiallelic"
    (candidate positions dropped because the VCF has records of more than one allele there)}.
    """
    cand = {str(c): np.unique(g["Pos"].to_numpy("int64")) for c, g in candidates.groupby("Chr", observed=True)}
    fmt, chunks = _chunks(path)
    stats = {"format": fmt, "records": 0, "sites": 0, "chroms": set(), "multiallelic": 0}
    found = []
    for chrom, pos, ref, alt in chunks:
        stats["records"] += len(chrom)
        if fmt == "vcf":
            snv = (ref.str.fullmatch(_SNV) & alt.str.fullmatch(_SNV)).to_numpy(bool)
        else:
            snv = np.ones(len(chrom), bool)
        stats["sites"] += int(snv.sum())
        positions = _positions(pos, path, snv)
        chroms = chrom.to_numpy(str)
        stats["chroms"].update(str(c) for c in pd.unique(chroms[snv]))
        # every record at a candidate position, the dropped ones too: a split multi-allelic site has a record of
        # another allele (an SNV, an indel or `*`) at the same position
        for c in pd.unique(chroms):
            if str(c) not in cand:
                continue
            rows = np.flatnonzero(chroms == c)
            rows = rows[np.isin(positions[rows], cand[str(c)])]
            if not len(rows):
                continue
            hit = pd.DataFrame({"Chr": str(c), "Pos": positions[rows], "SNV": snv[rows]})
            if fmt == "vcf":
                hit["Allele"] = (ref.iloc[rows].str.upper() + ">" + alt.iloc[rows].str.upper()).to_numpy()
            found.append(hit)
    if not found:
        return pd.DataFrame({"Chr": pd.Series(dtype=str), "Pos": pd.Series(dtype="int64")}), stats
    hits = pd.concat(found, ignore_index=True).drop_duplicates()  # an identical record listed twice counts once
    if fmt == "vcf":
        site = hits.groupby(["Chr", "Pos"], sort=False).agg(alleles=("Allele", "size"), snv=("SNV", "all"))
        multi = site["alleles"] > 1
        stats["multiallelic"] = int(multi.sum())
        out = site.index[(site["snv"] & ~multi).to_numpy(bool)].to_frame(index=False)
    else:
        out = hits[["Chr", "Pos"]]
    return out.astype({"Pos": "int64"}).reset_index(drop=True), stats


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
    multi = f" ({st['multiallelic']} multi-allelic candidate positions dropped)" if st["multiallelic"] else ""
    log.info("SNP sites %s (%s): %d records, %d %s, %d at a candidate position%s", path, st["format"].upper(),
             st["records"], st["sites"], what, len(matched), multi)
    if not st["sites"]:
        raise ValueError(f"--snp-sites {path} has no {what} ({st['records']} records): wrong file?")
    chroms = {str(c) for c in chroms}
    if not st["chroms"] & chroms:
        raise ValueError(f"--snp-sites {path} has no site on the chromosomes of the BED ({_shown(chroms)}; sites on "
                         f"{_shown(st['chroms'])}): chromosome naming or genome build of the sites file?")
    keys = pd.MultiIndex.from_frame(matched.astype({"Chr": str}))
    keep = pd.MultiIndex.from_frame(snp_df[["Chr", "Pos"]].astype({"Chr": str, "Pos": "int64"})).isin(keys)
    return snp_df.loc[keep].reset_index(drop=True)
