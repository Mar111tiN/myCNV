"""Step B/C: pileup (sample) or PONmatrix (PON) -> 10 bp coverage bins + het-SNP / pooled-SNP tables.

Both wrap the mawk chain `filterBed | pile2CNV` (sample) or `filterBed | PON2CNV` (PON) and then merge the GC ratio
and the mappability tracks. In the fixed coordinates (default) the bin ExonPos is taken from the BED
(`bed.bin_exon_pos`), so sample and PON bins get the same (Chr, Pos, ExonPos) whatever bases their streams hold.
"""

from __future__ import annotations

import os
import shlex
import tempfile

import pandas as pd

from .bed import bin_exon_pos, normalize_chrom, read_bed
from .tools import cat_cmd, cmd2df, log, mawk
from .tracks import add_gc, add_genmap


def _chrom_list(chrom) -> list[str] | None:
    """None (whole BED), one name or a list -> de-duplicated list in the given order (filterBed naming)."""
    if chrom is None:
        return None
    chroms = [chrom] if isinstance(chrom, str) else list(chrom)
    return list(dict.fromkeys(normalize_chrom(c) for c in chroms)) or None


def _filter_cmd(bed: str, chroms: list[str] | None, legacy_coords: bool) -> str:
    cmd = f"{mawk('filterBed')} {shlex.quote(bed)} -x"
    if chroms:
        cmd += f" -c {shlex.quote(','.join(chroms))}"
    if legacy_coords:
        cmd += " -L"
    return cmd


def _regions(bed: str, chroms: list[str] | None):
    regions = read_bed(bed, chroms)  # validates the BED order (filterBed would drop positions silently)
    if not regions:
        raise ValueError(f"BED {bed} has no region{' on ' + ', '.join(chroms) if chroms else ''}")
    return regions


def _fix_bins(cov_df: pd.DataFrame, regions, bin_size: int, legacy_coords: bool, source: str) -> pd.DataFrame:
    if cov_df.empty:
        raise ValueError(f"no coverage bins from {source}: the input is empty or header-only, or none of its positions "
                         "lies in the BED (chromosome naming, --chrom, or a BED of another kit)")
    if legacy_coords:
        return cov_df
    exon = bin_exon_pos(cov_df, regions, bin_size)
    streamed = cov_df["ExonPos"].to_numpy()
    if (exon > streamed).any():  # filterBed and bed.read_bed disagree: never expected
        raise RuntimeError("BED-derived bin ExonPos exceeds the ExonPos of a streamed base")
    cov_df = cov_df.copy()
    cov_df["ExonPos"] = exon
    return cov_df


def tn2cnv(
    pileup: str,
    bed: str,
    gc: str | None = None,
    genmap: str | None = None,
    chrom=None,
    normal_vaf=(0.2, 0.8),
    min_depth: int = 15,
    min_cov: float = 0,
    bin_size: int = 10,
    clean: bool = True,
    legacy_coords: bool = False,
    tmp_dir: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage bins and het SNPs of a paired pileup (normal = sample 1).

    pileup: `samtools mpileup ... NORMAL TUMOR | cleanpileup -d` output (clean=True, the somVarWES pileup files),
            or raw `samtools mpileup` text (clean=False; cleanpileup -d is added). Plain or gzipped.
    chrom:  restrict to one or more chromosomes (optional split, plan D16); default is the whole sample. GC
            normalisation needs the whole sample, so per-chromosome outputs must all go into one `combine`.
    legacy_coords: legacy bin labels and ExonPos (only to match PON coverage built by the legacy tools).
    Raises ValueError when no coverage bin is made (empty pileup, BED or --chrom do not match the pileup).
    """
    chroms = _chrom_list(chrom)
    regions = _regions(bed, chroms)
    read_cmd = cat_cmd(pileup)
    if not clean:
        read_cmd += f" | {mawk('cleanpileup')} -d"
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        snp_file = os.path.join(tmp, "sample.snp")
        cnv_cmd = (
            f"{mawk('pile2CNV')} -x -o {shlex.quote(snp_file)} -v {normal_vaf[0]} -V {normal_vaf[1]} "
            f"-d {min_depth} -c {min_cov} -b {bin_size}"
        )
        if legacy_coords:
            cnv_cmd += " -L"
        cmd = f"{read_cmd} | {_filter_cmd(bed, chroms, legacy_coords)} | {cnv_cmd}"
        log.info("Coverage and het SNPs from %s%s", pileup, f" ({', '.join(chroms)})" if chroms else "")
        cov_df = cmd2df(cmd)
        snp_df = pd.read_csv(snp_file, sep="\t")
    cov_df = _fix_bins(cov_df, regions, bin_size, legacy_coords, pileup)
    log.info("%d coverage bins, %d het SNPs", len(cov_df), len(snp_df))
    if snp_df.empty:
        log.warning("No het SNP in %s (normal VAF %s-%s, depth >= %s in every sample): no BAF track", pileup,
                    normal_vaf[0], normal_vaf[1], min_depth)
    cov_df = add_gc(cov_df, gc)
    cov_df, snp_df = add_genmap(cov_df, snp_df, track=genmap)
    return cov_df, snp_df


def pon2cnv(
    matrix: str,
    bed: str,
    gc: str | None = None,
    genmap: str | None = None,
    chrom=None,
    min_vaf: float = 0.05,
    min_depth: int = 50,
    min_cov: float = 0,
    bin_size: int = 10,
    legacy_coords: bool = False,
    tmp_dir: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample coverage bins and pooled SNP VAFs of one PONmatrix file (`matrix/{chrom}.pon.gz`)."""
    chroms = _chrom_list(chrom)
    regions = _regions(bed, chroms)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        snp_file = os.path.join(tmp, "pon.snp")
        cnv_cmd = (f"{mawk('PON2CNV')} -x -o {shlex.quote(snp_file)} -v {min_vaf} -d {min_depth} -c {min_cov} "
                   f"-b {bin_size}")
        if legacy_coords:
            cnv_cmd += " -L"
        cmd = f"{cat_cmd(matrix)} | {_filter_cmd(bed, chroms, legacy_coords)} | {cnv_cmd}"
        log.info("PON coverage and SNPs from %s", matrix)
        cov_df = cmd2df(cmd)
        snp_df = pd.read_csv(snp_file, sep="\t")
    cov_df = _fix_bins(cov_df, regions, bin_size, legacy_coords, matrix)
    cov_df = add_gc(cov_df, gc)
    # the CNV PON does not use the map columns: genmap is optional here
    cov_df, snp_df = add_genmap(cov_df, snp_df, track=genmap, required=False, names=("PON bins", "PON SNPs"))
    return cov_df, snp_df
