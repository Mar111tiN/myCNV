"""Step B/C: pileup (sample) or PONmatrix (PON) -> 10 bp coverage bins + het-SNP / pooled-SNP tables.

Both wrap the mawk chain `filterBed | pile2CNV` (sample) or `filterBed | PON2CNV` (PON) and then merge the GC ratio
and the mappability tracks. In the fixed coordinates (default) the bin ExonPos is taken from the BED
(`bed.bin_exon_pos`), so sample and PON bins get the same (Chr, Pos, ExonPos) whatever bases their streams hold.
A pair (`tn2cnv`, normal = sample 1) and a tumor-only sample (`tumor2cnv`, one-sample pileup, het-SNP candidates at
given common-SNP sites) run the same `filterBed | pile2CNV` pass.
"""

from __future__ import annotations

import gzip
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


def _regions(bed: str, chroms: list[str] | None, legacy_coords: bool = False):
    # validates the BED order (filterBed would drop positions silently); fixed coordinates merge overlapping/touching
    # regions exactly like filterBed without -L, so bin ExonPos and the streamed ExonPos agree
    regions = read_bed(bed, chroms, merge=not legacy_coords)
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


def _layout_error(pileup, clean: bool, got: str) -> ValueError:
    if clean:
        return ValueError(f"the pileup {pileup} does not start with the `cleanpileup -d` header (Chr Start Ref Depth1 "
                          f"Read1 ...): got {got}; a raw samtools mpileup needs --raw-pileup")
    return ValueError(f"the pileup {pileup} is not raw samtools mpileup text (3 + 3 per sample columns, got {got}); a "
                      "clean pileup (with header) is read without --raw-pileup")


def pileup_samples(pileup: str, clean: bool = True) -> int:
    """Number of samples of a pileup, from its first line.

    clean: the `cleanpileup -d` header `Chr Start Ref Depth1 Read1 [Depth2 Read2 ...]` (what `ebscore pileup` writes);
    raw `samtools mpileup` text (clean=False): 3 + 3 per sample columns. Raises ValueError for an empty file or a
    first line of another layout (e.g. a raw mpileup given without --raw-pileup).
    Reads (consumes) the first line: only for a regular file, not for a pipe or FIFO (`tn2cnv` checks those from the
    stream instead).
    """
    try:
        with (gzip.open(pileup, "rt") if str(pileup).endswith(".gz") else open(pileup)) as fh:
            first = fh.readline().rstrip("\r\n")
    except FileNotFoundError:
        raise
    except (OSError, EOFError, UnicodeDecodeError) as err:
        raise ValueError(f"cannot read the pileup {pileup}: {err}") from None
    if not first:
        raise ValueError(f"the pileup {pileup} is empty")
    fields = first.split("\t")
    if clean:
        n = (len(fields) - 3) // 2
        expected = ["Chr", "Start", "Ref"] + [f"{c}{k}" for k in range(1, n + 1) for c in ("Depth", "Read")]
        if n < 1 or fields != expected:
            raise _layout_error(pileup, clean, repr(first[:80]))
        return n
    if len(fields) < 6 or (len(fields) - 3) % 3 or fields[:2] == ["Chr", "Start"]:
        raise _layout_error(pileup, clean, str(len(fields)))
    return (len(fields) - 3) // 3


def _stream_samples(cov_df: pd.DataFrame, pileup, clean: bool) -> int:
    """Number of samples from the pile2CNV header `Chr Pos ExonPos Cov1 ... CovN`. A pileup without a valid header
    gives another header: filterBed passes its first line on as data, and pile2CNV takes that line as the header."""
    cols = [str(c) for c in cov_df.columns]
    n = len(cols) - 3
    if n < 1 or cols != ["Chr", "Pos", "ExonPos"] + [f"Cov{k}" for k in range(1, n + 1)]:
        raise _layout_error(pileup, clean, "the stream header " + repr("\t".join(cols)[:80]))
    return n


def _check_samples(pileup, n: int, tumor_only: bool) -> None:
    """A pair needs 2 samples (NORMAL TUMOR), a tumor-only sample 1; the error names the flag."""
    if tumor_only and n != 1:
        raise ValueError(f"--tumor-only needs a one-sample (TUMOR) pileup; {pileup} has {n} samples"
                         + (" (a NORMAL TUMOR pair: run without --tumor-only)" if n == 2 else ""))
    if not tumor_only and n != 2:
        raise ValueError(f"the pileup {pileup} has {n} sample{'s' if n != 1 else ''}; a pair needs NORMAL TUMOR (2)"
                         + (": a one-sample tumor pileup needs --tumor-only --snp-sites" if n == 1 else ""))


def _pile2cnv(pileup: str, bed: str, chroms, regions, vaf, min_depth: int, min_cov: float, bin_size: int, clean: bool,
              legacy_coords: bool, tmp_dir: str | None, tumor_only: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    """`[cleanpileup -d |] filterBed | pile2CNV` on a pileup: (coverage bins with the BED-derived ExonPos, het-SNP
    table). The het SNPs are the rows where sample 1 has an alt read, every sample has depth >= min_depth and the
    VAF of sample 1 is within `vaf` (inclusive): the normal of a pair, the tumor of a tumor-only sample.

    The sample count (2 for a pair, 1 for `tumor_only`) is checked on the first line of a regular file before the pass
    (an early error) and on the pile2CNV header after it. A pipe or FIFO (`/dev/stdin`, process substitution) cannot be
    read twice, so it gets only the second check: reading its first line here would take it away from the pass.
    """
    if os.path.isfile(pileup):
        _check_samples(pileup, pileup_samples(pileup, clean), tumor_only)
    elif not os.path.exists(pileup):
        raise FileNotFoundError(f"the pileup {pileup} does not exist")
    read_cmd = cat_cmd(pileup)
    if not clean:
        read_cmd += f" | {mawk('cleanpileup')} -d"
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        snp_file = os.path.join(tmp, "sample.snp")
        cnv_cmd = (
            f"{mawk('pile2CNV')} -x -o {shlex.quote(snp_file)} -v {vaf[0]} -V {vaf[1]} "
            f"-d {min_depth} -c {min_cov} -b {bin_size}"
        )
        if legacy_coords:
            cnv_cmd += " -L"
        cmd = f"{read_cmd} | {_filter_cmd(bed, chroms, legacy_coords)} | {cnv_cmd}"
        log.info("Coverage and het SNPs from %s%s", pileup, f" ({', '.join(chroms)})" if chroms else "")
        try:
            cov_df = cmd2df(cmd)
        except pd.errors.EmptyDataError:  # no line at all (a regular file is caught by pileup_samples)
            raise ValueError(f"the pileup {pileup} is empty") from None
        _check_samples(pileup, _stream_samples(cov_df, pileup, clean), tumor_only)
        snp_df = pd.read_csv(snp_file, sep="\t")
    cov_df = _fix_bins(cov_df, regions, bin_size, legacy_coords, pileup)
    return cov_df, snp_df


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
            or raw `samtools mpileup` text (clean=False; cleanpileup -d is added). Plain or gzipped. It must hold
            exactly 2 samples (a one-sample tumor pileup goes to `tumor2cnv`).
    chrom:  restrict to one or more chromosomes (optional split, plan D16); default is the whole sample. GC
            normalisation needs the whole sample, so per-chromosome outputs must all go into one `combine`.
    legacy_coords: legacy bin labels and ExonPos (only to match PON coverage built by the legacy tools).
    Raises ValueError when no coverage bin is made (empty pileup, BED or --chrom do not match the pileup).
    """
    chroms = _chrom_list(chrom)
    regions = _regions(bed, chroms, legacy_coords)
    cov_df, snp_df = _pile2cnv(pileup, bed, chroms, regions, normal_vaf, min_depth, min_cov, bin_size, clean,
                               legacy_coords, tmp_dir, tumor_only=False)
    log.info("%d coverage bins, %d het SNPs", len(cov_df), len(snp_df))
    if snp_df.empty:
        log.warning("No het SNP in %s (normal VAF %s-%s, depth >= %s in every sample): no BAF track", pileup,
                    normal_vaf[0], normal_vaf[1], min_depth)
    cov_df = add_gc(cov_df, gc)
    cov_df, snp_df = add_genmap(cov_df, snp_df, track=genmap)
    return cov_df, snp_df


def vaf_band(vaf, key: str = "hetSNP.tumorOnlyVAF") -> tuple[float, float]:
    """A VAF range [min, max] with 0 <= min <= max <= 1 (ValueError otherwise)."""
    try:
        lo, hi = (float(v) for v in vaf)
    except (TypeError, ValueError):
        raise ValueError(f"{key} must be two numbers [min, max], got {vaf!r}") from None
    if not 0 <= lo <= hi <= 1:
        raise ValueError(f"{key} needs 0 <= min <= max <= 1, got [{lo}, {hi}]")
    return lo, hi


def tumor2cnv(
    pileup: str,
    bed: str,
    snp_sites: str,
    gc: str | None = None,
    genmap: str | None = None,
    chrom=None,
    tumor_vaf=(0.1, 0.9),
    min_depth: int = 15,
    min_cov: float = 0,
    bin_size: int = 10,
    clean: bool = True,
    legacy_coords: bool = False,
    tmp_dir: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage bins and het-SNP candidates of a tumor-only sample (one-sample pileup, tumor = sample 1; plan D27).

    The coverage bins come from the same `filterBed | pile2CNV` pass as a pair's (Cov1 = tumor). Without a normal,
    the het-SNP candidates are the positions of `snp_sites` (VCF or 2-column TSV, `mycnv.sites`; N-Q3) where the tumor
    has depth >= min_depth and a VAF within `tumor_vaf` (config hetSNP.tumorOnlyVAF), with pile2CNV's alt reads
    (A/C/G/T, D/d, I/i). PON artefact sites (maxPONVAF) are dropped downstream, as for pairs.
    Returns (Chr Pos ExonPos GCratio map* Cov1, Chr Pos ExonPos map* Depth1 VAF1).
    """
    from .sites import check_sites, restrict_to_sites

    band = vaf_band(tumor_vaf)
    chroms = _chrom_list(chrom)
    regions = _regions(bed, chroms, legacy_coords)
    check_sites(snp_sites)  # before the pileup pass: missing, unreadable layout, no record
    cov_df, snp_df = _pile2cnv(pileup, bed, chroms, regions, band, min_depth, min_cov, bin_size, clean, legacy_coords,
                               tmp_dir, tumor_only=True)
    n_cand = len(snp_df)
    # the naming check uses all BED chromosomes: a --chrom split of a small panel may have no site on its chromosome
    bed_chroms = regions.keys() if chroms is None else read_bed(bed, merge=not legacy_coords).keys()
    snp_df = restrict_to_sites(snp_df, snp_sites, bed_chroms)
    log.info("%d coverage bins; tumor-only het SNPs: %d of %d tumor candidates (VAF %s-%s, depth >= %s) are SNP sites",
             len(cov_df), len(snp_df), n_cand, band[0], band[1], min_depth)
    if snp_df.empty:
        log.warning("No het SNP in %s at the sites of %s (tumor VAF %s-%s, depth >= %s): no BAF track", pileup,
                    snp_sites, band[0], band[1], min_depth)
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
    regions = _regions(bed, chroms, legacy_coords)
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
