"""`mycnv` command line: one subcommand per legacy step, plus `sample` for the whole per-sample chain."""

from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from .config import dump_config, load_config
from .tools import PRIMARY_CHROMS, log, setup_logging


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", help="YAML with parameter overrides (structure of `mycnv config`)")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging (prints the shell pipelines)")


def _tracks(p: argparse.ArgumentParser) -> None:
    p.add_argument("--bed", required=True, help="capture BED (the same flavour for samples and PON, plan L9)")
    p.add_argument("--gc", help="GC track: directory with {chrom}.gc100-10.gz or one file with a Chr column")
    p.add_argument("--genmap", help="mappability: directory with hg38_genmap.{chrom}.txt.gz or one file")
    p.add_argument("--legacy-coords", action="store_true",
                   help="legacy bin labels and ExonPos; only to combine with PON coverage built by the legacy tools")


def _sample_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--normal-vaf", nargs=2, type=float, metavar=("MIN", "MAX"), help="het-SNP normal VAF range (0.2 0.8)")
    p.add_argument("--min-depth", type=int, help="het-SNP minimum depth in every sample (15)")
    p.add_argument("--min-cov", type=float, help="minimum binned coverage for output (0)")
    p.add_argument("--raw-pileup", action="store_true", help="input is raw samtools mpileup text (cleanpileup -d is added)")
    p.add_argument("--chrom", help="restrict to one chromosome (optional split); default: whole sample")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mycnv", description=__doc__)
    parser.add_argument("--version", action="version", version=f"mycnv {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("config", help="print the effective parameters (defaults + --config)")
    _common(p)

    p = sub.add_parser("pon-coverage", help="PONmatrix -> PON cov/{chrom}.cov.gz + snp/{chrom}.snp.gz (per chromosome)")
    _common(p)
    _tracks(p)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--matrix-dir", help="directory with {chrom}.pon.gz (PONmatrix `matrix/`)")
    src.add_argument("--matrix", nargs="+", help="PONmatrix file(s) named {chrom}.pon.gz")
    p.add_argument("--out-dir", required=True, help="PON directory to write cov/ and snp/ into")
    p.add_argument("--chrom", nargs="+", help="chromosomes to process (default: all matrix files)")
    p.add_argument("--min-vaf", type=float, help="pooled VAF threshold for PON SNPs (0.05)")
    p.add_argument("--min-depth", type=int, help="pooled depth threshold for PON SNPs (50)")
    p.add_argument("--min-cov", type=float, help="minimum binned coverage for output (0)")
    p.add_argument("--threads", type=int, default=1, help="chromosomes processed in parallel")

    p = sub.add_parser("pon-combine", help="PON cov/ + snp/ -> CNV/pon.cov.{full,fullGC,filter,filterGC}.gz, CNV/pon.snp.gz")
    _common(p)
    p.add_argument("--pon-dir", required=True)
    p.add_argument("--out-dir", help="where CNV/ is written (default: --pon-dir)")
    p.add_argument("--std-factor", type=float, help="outlier threshold in PON std units (2)")
    p.add_argument("--chroms", nargs="+", default=PRIMARY_CHROMS, help="chromosomes to combine (chr1..chr22 chrX)")
    p.add_argument("--legacy", action="store_true", help="reproduce main: filter stats not recomputed after outlier removal")

    p = sub.add_parser("tn2cnv", help="paired pileup -> coverage bins + het SNPs (whole sample by default)")
    _common(p)
    _tracks(p)
    _sample_opts(p)
    p.add_argument("--pileup", required=True, help="NORMAL TUMOR pileup after `cleanpileup -d` (plain or .gz)")
    p.add_argument("--out-cov", required=True)
    p.add_argument("--out-snp", required=True)

    p = sub.add_parser("combine", help="tn2cnv outputs + PON -> raw/{s}.cov.gz + raw/{s}.snp (log2ratio vs PON)")
    _common(p)
    p.add_argument("--cov", nargs="+", required=True, help="tn2cnv coverage file(s), e.g. one per chromosome")
    p.add_argument("--snp", nargs="+", required=True, help="tn2cnv SNP file(s)")
    p.add_argument("--pon-dir", required=True)
    p.add_argument("--out-cov", required=True)
    p.add_argument("--out-snp", required=True)
    p.add_argument("--no-gc-normalize", action="store_true", help="mean-normalise instead of GC-normalise")
    _plot_opts(p)

    p = sub.add_parser("rolling", help="raw tables -> data/{s}.cnv.snp.gz + data/{s}.cnv.gz (filter + rolling windows)")
    _common(p)
    p.add_argument("--cov", required=True, help="raw/{s}.cov.gz")
    p.add_argument("--snp", required=True, help="raw/{s}.snp")
    p.add_argument("--out-snp", required=True, help="data/{s}.cnv.snp.gz")
    p.add_argument("--out-cnv", required=True, help="data/{s}.cnv.gz")
    _plot_opts(p)

    p = sub.add_parser("ascat-input", help="data/{s}.cnv.snp.gz -> ASCAT BAF/logR tables")
    _common(p)
    p.add_argument("--cnv", required=True, help="data/{s}.cnv.snp.gz")
    p.add_argument("--sample", required=True, help="sample ID used in file names and as column name")
    p.add_argument("--tumor", help="tumor name (default: from <patient>_<tumor>-<normal>)")
    p.add_argument("--normal", help="normal name (default: from <patient>_<tumor>-<normal>)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--no-mock-germline", action="store_true")

    p = sub.add_parser("sample", help="whole per-sample chain: tn2cnv -> combine -> rolling -> ascat-input")
    _common(p)
    _tracks(p)
    _sample_opts(p)
    p.add_argument("--pileup", required=True)
    p.add_argument("--pon-dir", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--tumor")
    p.add_argument("--normal")
    p.add_argument("--out-dir", required=True, help="writes raw/, data/, pre/ (and plots/) here")
    p.add_argument("--no-gc-normalize", action="store_true")
    p.add_argument("--no-mock-germline", action="store_true")
    p.add_argument("--plots", action="store_true", help="write the legacy JPG plots into plots/")
    return parser


def _plot_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--plots-dir", help="write the legacy JPG plots here")
    p.add_argument("--sample", help="sample ID for plot names (with --plots-dir)")
    p.add_argument("--tumor")
    p.add_argument("--normal")


def _config(args) -> dict:
    over: dict = {}

    def put(path, value):
        if value is None:
            return
        d = over
        for key in path[:-1]:
            d = d.setdefault(key, {})
        d[path[-1]] = value

    put(("hetSNP", "normalVAF"), getattr(args, "normal_vaf", None))
    if args.cmd in ("tn2cnv", "sample"):
        put(("hetSNP", "minDepth"), getattr(args, "min_depth", None))
        put(("coverage", "minCov"), getattr(args, "min_cov", None))
    if args.cmd == "pon-coverage":
        put(("PONcoverage", "minVAF"), args.min_vaf)
        put(("PONcoverage", "minDepth"), args.min_depth)
        put(("PONcoverage", "minCov"), args.min_cov)
    if args.cmd == "pon-combine":
        put(("PONcoverage", "stdFactor"), args.std_factor)
    if getattr(args, "no_gc_normalize", False):
        put(("coverage", "GCnormalize"), False)
    return load_config(args.config, over)


def _names(args):
    from .convert import split_sample_name

    if args.tumor and args.normal:
        return args.tumor, args.normal
    return split_sample_name(args.sample)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    cfg = _config(args)

    from . import pipeline as pl
    from .tools import read_tsv, write_tsv

    if args.cmd == "config":
        sys.stdout.write(dump_config(cfg))
        return 0

    if args.cmd == "pon-coverage":
        if args.matrix_dir:
            files = pl.find_matrix_files(args.matrix_dir, args.chrom)
        else:
            files = {}
            for path in args.matrix:
                chrom = os.path.basename(path).split(".pon")[0]
                if not args.chrom or chrom in args.chrom:
                    files[chrom] = path
        if not files:
            log.error("no PONmatrix files found")
            return 1
        c = cfg["PONcoverage"]
        pl.pon_coverage(
            files, args.out_dir, threads=args.threads, verbose=args.verbose, bed=args.bed, gc=args.gc,
            genmap=args.genmap, min_vaf=c["minVAF"], min_depth=c["minDepth"], min_cov=c["minCov"],
            legacy_coords=args.legacy_coords,
        )
        return 0

    if args.cmd == "pon-combine":
        pl.pon_combine(args.pon_dir, args.out_dir, std_factor=cfg["PONcoverage"]["stdFactor"], chroms=args.chroms,
                       legacy=args.legacy)
        return 0

    if args.cmd == "tn2cnv":
        from .raw import tn2cnv

        sc, cc = cfg["hetSNP"], cfg["coverage"]
        cov_df, snp_df = tn2cnv(
            args.pileup, args.bed, gc=args.gc, genmap=args.genmap, chrom=args.chrom, normal_vaf=sc["normalVAF"],
            min_depth=sc["minDepth"], min_cov=cc["minCov"], clean=not args.raw_pileup, legacy_coords=args.legacy_coords,
        )
        write_tsv(cov_df, args.out_cov)
        write_tsv(snp_df, args.out_snp)
        return 0

    if args.cmd == "combine":
        from .combine import concat_tables

        cov_df, snp_df = pl.combine_step(concat_tables(args.cov), concat_tables(args.snp), args.pon_dir, cfg)
        write_tsv(cov_df, args.out_cov)
        write_tsv(snp_df, args.out_snp)
        if args.plots_dir:
            pl.raw_plots(cov_df, snp_df, args.sample, *_names(args), args.plots_dir)
        return 0

    if args.cmd == "rolling":
        roll_snp, cnv_df, _ = pl.rolling_step(read_tsv(args.cov), read_tsv(args.snp), cfg)
        write_tsv(roll_snp, args.out_snp)
        write_tsv(cnv_df, args.out_cnv)
        if args.plots_dir:
            pl.cnv_plots(roll_snp, cnv_df, args.sample, *_names(args), args.plots_dir)
        return 0

    if args.cmd == "ascat-input":
        from .convert import write_ascat

        write_ascat(read_tsv(args.cnv), args.sample, args.out_dir, tumor=args.tumor, normal=args.normal,
                    mock_germline=not args.no_mock_germline)
        return 0

    if args.cmd == "sample":
        pl.run_sample(
            args.pileup, args.bed, args.pon_dir, args.sample, args.out_dir, cfg, gc=args.gc, genmap=args.genmap,
            tumor=args.tumor, normal=args.normal, chrom=args.chrom, clean=not args.raw_pileup,
            legacy_coords=args.legacy_coords, plots=args.plots, mock_germline=not args.no_mock_germline,
        )
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
