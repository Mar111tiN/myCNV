"""`mycnv` command line: one subcommand per legacy step, plus `sample` for the whole per-sample chain.

Option vocabulary shared with ebscore: --ponmatrix DIR (the PONmatrix folder: matrix/ from ebscore, cov/ snp/ CNV/ from
mycnv), -o/--outdir DIR, repeatable --chrom, --bed, --pileup, -t/--threads, --tmpdir, -v/--verbose, --quiet. Inputs
are never written to: every command that writes a folder needs -o/--outdir.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

from . import __version__
from .config import dump_config, load_config
from .tools import PRIMARY_CHROMS, log, setup_logging, unique


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", help="YAML with parameter overrides (structure of `mycnv config`)")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="parallel processes (pon-coverage: chromosomes in parallel; elsewhere accepted, serial)")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging (prints the shell pipelines)")
    p.add_argument("--quiet", action="store_true", help="warnings and errors only")


def _chrom(p: argparse.ArgumentParser, help: str) -> None:
    p.add_argument("--chrom", action="append", metavar="CHROM", help=f"{help} (repeatable)")


def _tracks(p: argparse.ArgumentParser, genmap_required: bool = True) -> None:
    p.add_argument("--bed", required=True, help="capture BED: the kit's padded BED, the same file for samples and PON (plan D20)")
    p.add_argument("--gc", required=True, help="GC track: directory with {chrom}.gc100-10.gz or one file")
    p.add_argument("--genmap", required=genmap_required,
                   help="mappability: directory with hg38_genmap.{chrom}.txt.gz or one file"
                   + ("" if genmap_required else " (optional: the CNV PON does not use it)"))
    p.add_argument("--legacy-coords", action="store_true",
                   help="legacy bin labels and ExonPos; only to combine with PON coverage built by the legacy tools")
    p.add_argument("--tmpdir", help="scratch folder (default: system temp)")


def _sample_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--normal-vaf", nargs=2, type=float, metavar=("MIN", "MAX"), help="het-SNP normal VAF range (0.2 0.8)")
    p.add_argument("--min-depth", type=int, help="het-SNP minimum depth in every sample (15)")
    p.add_argument("--min-cov", type=float, help="minimum binned coverage for output (0)")
    p.add_argument("--raw-pileup", action="store_true", help="input is raw samtools mpileup text (cleanpileup -d is added)")


def _pon_match(p: argparse.ArgumentParser) -> None:
    p.add_argument("--min-pon-match", type=float, default=0.9, metavar="FRACTION",
                   help="fail if fewer sample bins have a PON bin (default 0.9)")


def _names(p: argparse.ArgumentParser) -> None:
    p.add_argument("--tumor-name", help="tumor label in file names (default: from <patient>_<tumor>-<normal>)")
    p.add_argument("--normal-name", help="normal label in file names (default: from <patient>_<tumor>-<normal>)")


def _plot_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--plots-dir", help="write the legacy JPG plots here (needs --sample)")
    p.add_argument("--sample", help="sample ID for plot names")
    _names(p)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mycnv", description=__doc__, allow_abbrev=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"mycnv {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True, metavar="COMMAND")

    def add(name, help):
        return sub.add_parser(name, help=help, description=help, allow_abbrev=False)

    p = add("config", "print the effective parameters (defaults + --config)")
    _common(p)

    p = add("tracks", "build the GC (gc100-10) or mappability (genmap) tracks from a FASTA (genome-bundle builds, D22)")
    tsub = p.add_subparsers(dest="track", required=True, metavar="TRACK")
    for name, help in [
        ("gc", "{chrom}.gc100-10.gz: GC/(GC+AT) over [S-49, S+50] on a 10 bp grid, 2 decimals"),
        ("genmap", "hg38_genmap.{chrom}.txt.gz: genmap mappability (1 / occurrences) per k_e mode, dense over --bed"),
    ]:
        t = tsub.add_parser(name, help=help, description=help, allow_abbrev=False)
        t.add_argument("--fasta", required=True, help="FASTA with .fai (plain or bgzipped)")
        t.add_argument("-o", "--outdir", required=True, help="output folder")
        t.add_argument("--chrom", action="append", metavar="CHROM",
                       help="sequences to write, as named in the FASTA (repeatable; default: all)")
        if name == "genmap":
            t.add_argument("--bed", help="selection: rows are written for these regions only (e.g. the kit's padded "
                                         "BED, or exons +-100 bp); default: whole chromosomes (one row per base)")
            t.add_argument("--modes", default="30_0,50_0,75_1,100_2",
                           help="comma-separated <k-mer length>_<errors> (default: the four mycnv reads)")
            t.add_argument("--index", help="genmap index of the WHOLE genome: used if it exists, else built there and "
                                           "kept (~6 bytes per base of RAM; genome-wide only as a kumo job or alone on a "
                                           ">= 64 GB host). Default: built in --tmpdir and removed")
            t.add_argument("--tmpdir", help="scratch folder (default: system temp)")
        _common(t)

    p = add("pon-coverage", "PONmatrix matrix/ -> PON cov/{chrom}.cov.gz + snp/{chrom}.snp.gz (per chromosome)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--ponmatrix", help="PONmatrix folder (reads matrix/{chrom}.pon.gz)")
    src.add_argument("--matrix", nargs="+", action="extend", help="PONmatrix file(s) named {chrom}.pon.gz (repeatable)")
    p.add_argument("-o", "--outdir", required=True, help="folder that gets cov/ and snp/")
    _chrom(p, "chromosomes to process (default: all matrix files)")
    _tracks(p, genmap_required=False)
    p.add_argument("--min-vaf", type=float, help="pooled VAF threshold for PON SNPs (0.05)")
    p.add_argument("--min-depth", type=int, help="pooled depth threshold for PON SNPs (50)")
    p.add_argument("--min-cov", type=float, help="minimum binned coverage for output (0)")
    _common(p)

    p = add("pon-combine", "PON cov/ + snp/ -> CNV/pon.cov.{full,fullGC,filter,filterGC}.gz, CNV/pon.snp.gz")
    p.add_argument("--ponmatrix", required=True, help="folder with cov/ and snp/ (the pon-coverage -o folder)")
    p.add_argument("-o", "--outdir", required=True, help="folder that gets CNV/ (may be the --ponmatrix folder)")
    _chrom(p, "chromosomes to combine (default: chr1..chr22, chrX)")
    p.add_argument("--std-factor", type=float, help="outlier threshold in PON std units (2)")
    p.add_argument("--legacy-stats", action="store_true",
                   help="reproduce main: filter stats not recomputed after outlier removal")
    _common(p)

    p = add("tn2cnv", "paired pileup -> coverage bins + het SNPs (whole sample by default)")
    p.add_argument("--pileup", required=True, help="NORMAL TUMOR pileup after `cleanpileup -d` (plain or .gz)")
    _tracks(p)
    _sample_opts(p)
    _chrom(p, "only these chromosomes (optional split; all outputs then go into one `combine`)")
    p.add_argument("--out-cov", required=True, help="coverage bins (.gz for gzip)")
    p.add_argument("--out-snp", required=True, help="het SNPs")
    _common(p)

    p = add("combine", "tn2cnv outputs + CNV PON -> raw/{s}.cov.gz + raw/{s}.snp (log2ratio vs PON)")
    p.add_argument("--cov", nargs="+", action="extend", required=True,
                   help="tn2cnv coverage file(s): the whole sample, e.g. one per chromosome (repeatable)")
    p.add_argument("--snp", nargs="+", action="extend", required=True, help="tn2cnv SNP file(s) (repeatable)")
    p.add_argument("--ponmatrix", required=True, help="CNV PON folder (reads CNV/)")
    p.add_argument("--out-cov", required=True)
    p.add_argument("--out-snp", required=True)
    p.add_argument("--no-gc-normalize", action="store_true", help="mean-normalise instead of GC-normalise")
    p.add_argument("--legacy-coords", action="store_true", help="the tn2cnv outputs were made with --legacy-coords")
    _pon_match(p)
    p.add_argument("--allow-missing-chroms", action="store_true",
                   help="do not fail when PON chromosomes have no sample bins (partial sample: wrong normalisation)")
    _plot_opts(p)
    _common(p)

    p = add("rolling", "raw tables -> data/{s}.cnv.snp.gz + data/{s}.cnv.gz (filter + rolling windows)")
    p.add_argument("--cov", required=True, help="raw/{s}.cov.gz")
    p.add_argument("--snp", required=True, help="raw/{s}.snp")
    p.add_argument("--out-snp", required=True, help="data/{s}.cnv.snp.gz")
    p.add_argument("--out-cnv", required=True, help="data/{s}.cnv.gz")
    _plot_opts(p)
    _common(p)

    p = add("ascat-input", "data/{s}.cnv.snp.gz -> ASCAT BAF/logR tables")
    p.add_argument("--cnv", required=True, help="data/{s}.cnv.snp.gz")
    p.add_argument("--sample", required=True, help="sample ID used in file names and as column name")
    _names(p)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--no-mock-germline", action="store_true")
    p.add_argument("--allow-empty", action="store_true", help="write the tables even if there is no SNP or logR")
    _common(p)

    p = add("sample", "whole per-sample chain: tn2cnv -> combine -> rolling -> ascat-input (always the whole sample)")
    p.add_argument("--pileup", required=True, help="NORMAL TUMOR pileup after `cleanpileup -d` (plain or .gz)")
    p.add_argument("--ponmatrix", required=True, help="CNV PON folder (reads CNV/)")
    p.add_argument("--sample", required=True, help="sample ID used in file names")
    _names(p)
    p.add_argument("-o", "--outdir", required=True, help="writes raw/, data/, pre/ (and plots/) here")
    _tracks(p)
    _sample_opts(p)
    p.add_argument("--no-gc-normalize", action="store_true", help="mean-normalise instead of GC-normalise")
    p.add_argument("--no-mock-germline", action="store_true")
    _pon_match(p)
    p.add_argument("--plots", action="store_true", help="write the legacy JPG plots into plots/")
    p.add_argument("--allow-empty", action="store_true",
                   help="write header-only ASCAT tables when no SNP passes the filter (default: fail)")
    _common(p)
    return parser


def _config(args) -> dict:
    over: dict = {}

    def put(path, value):
        if value is None:
            return
        d = over
        for key in path[:-1]:
            d = d.setdefault(key, {})
        d[path[-1]] = value

    if args.cmd in ("tn2cnv", "sample"):
        put(("hetSNP", "normalVAF"), args.normal_vaf)
        put(("hetSNP", "minDepth"), args.min_depth)
        put(("coverage", "minCov"), args.min_cov)
    if args.cmd == "pon-coverage":
        put(("PONcoverage", "minVAF"), args.min_vaf)
        put(("PONcoverage", "minDepth"), args.min_depth)
        put(("PONcoverage", "minCov"), args.min_cov)
    if args.cmd == "pon-combine":
        put(("PONcoverage", "stdFactor"), args.std_factor)
    if getattr(args, "no_gc_normalize", False):
        put(("coverage", "GCnormalize"), False)
    return load_config(args.config, over)


def _check(parser: argparse.ArgumentParser, args) -> None:
    if args.cmd == "tracks":  # FASTA sequence names are used as given (no chr prefix rule)
        if args.chrom:
            args.chrom = unique(args.chrom)
    elif getattr(args, "chrom", None):
        from .bed import normalize_chrom

        args.chrom = unique(normalize_chrom(c) for c in args.chrom)
    if getattr(args, "plots_dir", None) and not args.sample:
        parser.error("--plots-dir needs --sample")
    if (getattr(args, "tumor_name", None) is None) != (getattr(args, "normal_name", None) is None):
        parser.error("give both --tumor-name and --normal-name, or neither")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if hasattr(args, "min_pon_match") and not 0 <= args.min_pon_match <= 1:
        parser.error("--min-pon-match must be within 0..1")


def _plot_names(args):
    from .convert import resolve_names

    return resolve_names(args.sample, args.tumor_name, args.normal_name)


def _run(args, cfg) -> None:
    from . import pipeline as pl
    from .tools import read_tsv, write_tsv

    if args.cmd == "config":
        sys.stdout.write(dump_config(cfg))
    elif args.cmd == "tracks":
        from . import trackgen

        if args.track == "gc":
            trackgen.make_gc_track(args.fasta, args.outdir, chroms=args.chrom, threads=args.threads)
        else:
            index = args.index
            if index and not os.path.exists(index):
                trackgen.build_genmap_index(args.fasta, index)
            trackgen.make_genmap_track(args.fasta, args.outdir, chroms=args.chrom, bed=args.bed, modes=args.modes,
                                       index=index, threads=args.threads, tmp_dir=args.tmpdir)
    elif args.cmd == "pon-coverage":
        if args.ponmatrix:
            files = pl.find_matrix_files(os.path.join(args.ponmatrix, "matrix"), args.chrom)
        else:
            files = pl.matrix_files_from_paths(args.matrix, args.chrom)
        if not files:
            raise FileNotFoundError("no PONmatrix files {chrom}.pon.gz found"
                                    + (f" for {', '.join(args.chrom)}" if args.chrom else ""))
        c = cfg["PONcoverage"]
        pl.pon_coverage(
            files, args.outdir, threads=args.threads, verbose=args.verbose, quiet=args.quiet, bed=args.bed,
            gc=args.gc, genmap=args.genmap, min_vaf=c["minVAF"], min_depth=c["minDepth"], min_cov=c["minCov"],
            legacy_coords=args.legacy_coords, tmp_dir=args.tmpdir,
        )
    elif args.cmd == "pon-combine":
        pl.pon_combine(args.ponmatrix, args.outdir, std_factor=cfg["PONcoverage"]["stdFactor"],
                       chroms=args.chrom or PRIMARY_CHROMS, legacy=args.legacy_stats)
    elif args.cmd == "tn2cnv":
        from .raw import tn2cnv

        sc, cc = cfg["hetSNP"], cfg["coverage"]
        cov_df, snp_df = tn2cnv(
            args.pileup, args.bed, gc=args.gc, genmap=args.genmap, chrom=args.chrom, normal_vaf=sc["normalVAF"],
            min_depth=sc["minDepth"], min_cov=cc["minCov"], clean=not args.raw_pileup, legacy_coords=args.legacy_coords,
            tmp_dir=args.tmpdir,
        )
        write_tsv(cov_df, args.out_cov)
        write_tsv(snp_df, args.out_snp)
    elif args.cmd == "combine":
        from .combine import concat_tables

        names = _plot_names(args) if args.plots_dir else None
        cov_df, snp_df = pl.combine_step(concat_tables(args.cov), concat_tables(args.snp), args.ponmatrix, cfg,
                                         legacy_coords=args.legacy_coords, min_pon_match=args.min_pon_match,
                                         allow_missing_chroms=args.allow_missing_chroms)
        write_tsv(cov_df, args.out_cov)
        write_tsv(snp_df, args.out_snp)
        if args.plots_dir:
            pl.raw_plots(cov_df, snp_df, args.sample, *names, args.plots_dir)
    elif args.cmd == "rolling":
        names = _plot_names(args) if args.plots_dir else None
        roll_snp, cnv_df, _ = pl.rolling_step(read_tsv(args.cov), read_tsv(args.snp), cfg)
        write_tsv(roll_snp, args.out_snp)
        write_tsv(cnv_df, args.out_cnv)
        if args.plots_dir:
            pl.cnv_plots(roll_snp, cnv_df, args.sample, *names, args.plots_dir)
    elif args.cmd == "ascat-input":
        from .convert import write_ascat

        write_ascat(read_tsv(args.cnv), args.sample, args.outdir, tumor=args.tumor_name, normal=args.normal_name,
                    mock_germline=not args.no_mock_germline, allow_empty=args.allow_empty)
    elif args.cmd == "sample":
        pl.run_sample(
            args.pileup, args.bed, args.ponmatrix, args.sample, args.outdir, cfg, gc=args.gc, genmap=args.genmap,
            tumor=args.tumor_name, normal=args.normal_name, clean=not args.raw_pileup,
            legacy_coords=args.legacy_coords, plots=args.plots, mock_germline=not args.no_mock_germline,
            min_pon_match=args.min_pon_match, allow_empty=args.allow_empty, tmp_dir=args.tmpdir,
        )


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _check(parser, args)
    setup_logging(args.verbose, args.quiet)
    try:
        cfg = _config(args)
        _run(args, cfg)
    except (ValueError, FileNotFoundError) as err:
        log.error("%s", err)
        return 2
    except subprocess.CalledProcessError as err:
        log.error("shell pipeline failed with exit status %s (see the tool messages above): %s", err.returncode,
                  err.cmd[-1] if isinstance(err.cmd, list) else err.cmd)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
