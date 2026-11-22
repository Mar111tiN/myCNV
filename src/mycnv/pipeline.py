"""High-level steps used by the CLI (and usable from notebooks). Output layout follows somVarWES:

PON dir:    cov/{chrom}.cov.gz  snp/{chrom}.snp.gz  CNV/pon.cov.{full,fullGC,filter,filterGC}.gz  CNV/pon.snp.gz
sample dir: raw/{s}.cov.gz  raw/{s}.snp  data/{s}.cnv.snp.gz  data/{s}.cnv.gz  pre/{s}_*_{baf,logr}.tsv  plots/*.jpg
"""

from __future__ import annotations

import os
import re
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from . import combine as cmb
from .convert import split_sample_name, write_ascat
from .raw import pon2cnv, tn2cnv
from .rolling_cov import rolling_coverage
from .rolling_snp import remerge_cnv, rolling_snp
from .tools import PRIMARY_CHROMS, chrom_key, log, read_tsv, setup_logging, write_tsv


# ##################### PON ###########################################################
def find_matrix_files(matrix_dir: str, chroms=None) -> dict[str, str]:
    """{chrom: path} for matrix/{chrom}.pon.gz (or .pon) files, in natural chromosome order."""
    found = {}
    for name in os.listdir(matrix_dir):
        m = re.match(r"^(.+)\.pon(\.gz)?$", name)
        if m and (chroms is None or m.group(1) in chroms):
            found[m.group(1)] = os.path.join(matrix_dir, name)
    return dict(sorted(found.items(), key=lambda kv: chrom_key(kv[0])))


def _pon_coverage_one(args) -> tuple[str, int, int]:
    chrom, matrix, out_dir, kw, verbose = args
    setup_logging(verbose)
    cov_df, snp_df = pon2cnv(matrix, chrom=chrom, **kw)
    write_tsv(cov_df, os.path.join(out_dir, "cov", f"{chrom}.cov.gz"), compression="gzip")
    write_tsv(snp_df, os.path.join(out_dir, "snp", f"{chrom}.snp.gz"), compression="gzip")
    return chrom, len(cov_df), len(snp_df)


def pon_coverage(matrix_files: dict[str, str], out_dir: str, threads: int = 1, verbose: bool = False, **kw):
    """PON2CNV for every chromosome (one job per chromosome, `threads` in parallel)."""
    jobs = [(chrom, path, out_dir, kw, verbose) for chrom, path in matrix_files.items()]
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            results = list(pool.map(_pon_coverage_one, jobs))
    else:
        results = [_pon_coverage_one(j) for j in jobs]
    for chrom, n_cov, n_snp in results:
        log.info("%s: %d coverage bins, %d SNPs", chrom, n_cov, n_snp)
    return results


def pon_combine(pon_dir: str, out_dir: str | None = None, std_factor: float = 2, chroms=PRIMARY_CHROMS,
                legacy: bool = False) -> list[str]:
    tables = cmb.make_pon_coverage(pon_dir, std_factor=std_factor, chroms=chroms, legacy=legacy)
    snp_df = cmb.make_pon_snp(pon_dir, chroms=chroms)
    return cmb.write_pon(out_dir or pon_dir, tables, snp_df)


# ##################### SAMPLE ########################################################
def raw_plots(cov_df, snp_df, sample: str, tumor: str, normal: str, plot_dir: str) -> list[str]:
    from .plot import plot_cov, plot_snp
    from .plot_params import covfig_params, log1, log2, snpfig_params, vaf1, vaf2

    os.makedirs(plot_dir, exist_ok=True)
    base = os.path.join(plot_dir, sample)
    files = []
    for name, plots in [(tumor, [vaf2]), (normal, [vaf1])] if not snp_df.empty else []:
        fig, *_ = plot_snp(snp_df, plots=plots, chroms="all", region="", **snpfig_params)
        files.append(f"{base}_{name}.snp.raw.jpg")
        fig.savefig(files[-1])
    for name, plots in [(tumor, [log2]), (normal, [log1])] if not cov_df.empty else []:
        fig, *_ = plot_cov(cov_df, plots=plots, chroms="all", region="", **covfig_params)
        files.append(f"{base}_{name}.cov.raw.jpg")
        fig.savefig(files[-1])
    _close_all()
    return files


def cnv_plots(snp_df, cnv_df, sample: str, tumor: str, normal: str, plot_dir: str) -> list[str]:
    from .plot import plot_CNV
    from .plot_params import L2Rmean1, L2Rmean2, cnvfig_params, log1, log2, vaf1, vaf2

    os.makedirs(plot_dir, exist_ok=True)
    base = os.path.join(plot_dir, sample)
    files = []
    for df, suffix in [(cnv_df, "cnv"), (snp_df, "snp.cnv")]:
        if df.empty:
            log.warning("No rows for the %s plots of %s (e.g. no SNP passed the filter): skipped", suffix, sample)
            continue
        for name, cov_plots, snp_plots in [(tumor, [log2, L2Rmean2], [vaf2]), (normal, [log1, L2Rmean1], [vaf1])]:
            fig, *_ = plot_CNV(df, cov_plots=cov_plots, snp_plots=snp_plots, chroms="all", region="", **cnvfig_params)
            files.append(f"{base}_{name}.{suffix}.jpg")
            fig.savefig(files[-1])
    _close_all()
    return files


def _close_all():
    import matplotlib.pyplot as plt

    plt.close("all")


def combine_step(cov_df, snp_df, pon_dir: str, config: dict):
    return cmb.combine_sample(cov_df, snp_df, pon_dir, gc_normalize=config["coverage"]["GCnormalize"])


def rolling_step(cov_df: pd.DataFrame, snp_df: pd.DataFrame, config: dict):
    """raw tables -> (rolled SNP table = data/*.cnv.snp.gz, full CNV track = data/*.cnv.gz, rolled coverage)."""
    cov_df = cmb.filter_cov(cov_df, config)
    roll_cov = rolling_coverage(cov_df, config)
    roll_snp = rolling_snp(snp_df, roll_cov, config)
    cnv_df = remerge_cnv(roll_snp, roll_cov)
    return roll_snp, cnv_df, roll_cov


def sample_paths(out_dir: str, sample: str) -> dict[str, str]:
    return {
        "raw_cov": os.path.join(out_dir, "raw", f"{sample}.cov.gz"),
        "raw_snp": os.path.join(out_dir, "raw", f"{sample}.snp"),
        "cnv_snp": os.path.join(out_dir, "data", f"{sample}.cnv.snp.gz"),
        "cnv": os.path.join(out_dir, "data", f"{sample}.cnv.gz"),
        "pre": os.path.join(out_dir, "pre"),
        "plots": os.path.join(out_dir, "plots"),
    }


def run_sample(pileup: str, bed: str, pon_dir: str, sample: str, out_dir: str, config: dict, gc=None, genmap=None,
               tumor=None, normal=None, chrom=None, clean=True, legacy_coords=False, plots=False,
               mock_germline=True) -> dict[str, str]:
    """Whole per-sample chain: tn2cnv -> combine -> rolling -> ASCAT input (one task per sample, plan D16)."""
    if tumor is None or normal is None:
        tumor, normal = split_sample_name(sample)
    p = sample_paths(out_dir, sample)
    sc, cc = config["hetSNP"], config["coverage"]
    cov_df, snp_df = tn2cnv(
        pileup, bed, gc=gc, genmap=genmap, chrom=chrom, normal_vaf=sc["normalVAF"], min_depth=sc["minDepth"],
        min_cov=cc["minCov"], clean=clean, legacy_coords=legacy_coords,
    )
    cov_df, snp_df = combine_step(cov_df, snp_df, pon_dir, config)
    write_tsv(cov_df, p["raw_cov"], compression="gzip")
    write_tsv(snp_df, p["raw_snp"])
    if plots:
        raw_plots(cov_df, snp_df, sample, tumor, normal, p["plots"])
    roll_snp, cnv_df, _ = rolling_step(cov_df, snp_df, config)
    write_tsv(roll_snp, p["cnv_snp"], compression="gzip")
    write_tsv(cnv_df, p["cnv"], compression="gzip")
    if plots:
        cnv_plots(roll_snp, cnv_df, sample, tumor, normal, p["plots"])
    write_ascat(roll_snp, sample, p["pre"], tumor=tumor, normal=normal, mock_germline=mock_germline)
    return p


def read_raw(cov: str, snp: str):
    return read_tsv(cov), read_tsv(snp)
