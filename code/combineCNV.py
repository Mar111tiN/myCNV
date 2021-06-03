# try with pd merge
import os
import numpy as np
import pandas as pd
from script_utils_CNV import show_output


def make_PON_snp(config={}, save=True):
    """
    combine the PON coverages for all chroms
    save to hardcoded place in PON_path
    """

    # paths
    pon_path = config["PON_path"]

    chrom_list = [f"chr{c + 1}" for c in range(22)] + ["chrX"]
    dfs = []
    for chrom in chrom_list:
        PON_snp_file = os.path.join(pon_path, f"snp/{chrom}.snp.gz")
        if os.path.isfile(PON_snp_file):
            show_output(f"Loading PON SNP file {PON_snp_file}")
            snp_df = pd.read_csv(PON_snp_file, sep="\t", compression="gzip")
            dfs.append(snp_df)
        else:
            show_output(
                f"Could not find PON coverage file {PON_snp_file}", color="warning"
            )
    snp_df = (
        pd.concat(dfs)
        .reset_index(drop=True)
        .rename({"VAF": "PONVAF", "Depth": "PONDepth"}, axis=1)
    )
    # make chrom categorical
    snp_df.loc[:, "Chr"] = pd.Categorical(snp_df["Chr"], chrom_list)

    # save file
    if save:
        PON_snp_file = os.path.join(pon_path, "CNV/pon.snp.gz")
        show_output(f"Saving combined PON SNP file {PON_snp_file}.")
        snp_df.to_csv(PON_snp_file, sep="\t", index=False, compression="gzip")
        show_output("Finished", color="success")
    return snp_df


def combine_PON_coverage(config={}):
    """
    combine the PON coverages for all chroms
    """

    # paths
    pon_path = config["PON_path"]

    chrom_list = [f"chr{c + 1}" for c in range(22)] + ["chrX"]
    dfs = []
    for chrom in chrom_list:
        PON_cov_file = os.path.join(pon_path, f"cov/{chrom}.cov.gz")
        if os.path.isfile(PON_cov_file):
            show_output(f"Loading coverage file {PON_cov_file}")
            cov_df = pd.read_csv(PON_cov_file, sep="\t", compression="gzip")
            dfs.append(cov_df)
        else:
            show_output(f"Could not find PON coverage file {PON_cov_file}")
    cov_df = pd.concat(dfs).reset_index(drop=True)
    # make chrom categorical
    cov_df.loc[:, "Chr"] = pd.Categorical(cov_df["Chr"], chrom_list)
    return cov_df


def normalize_GC_col(cov_df, col):
    """
    normalizes one coverage column for GC ratio
    """
    # compute the normalizer df
    # for each GCratio, norm_df has the difference of the respective mean from arbitrary norm coverage 100
    # remove chrX for the normalization or male genomes will have slightly greater mean
    norm_df = (
        (100 / cov_df.query('Chr != "chrX"').groupby("GCratio").agg({col: "mean"}))
        .reset_index()
        .rename({col: "factor"}, axis=1)
    )

    # merge to get the factor
    cov_df = cov_df.merge(norm_df)
    # adjust coverage using the factor from norm_df
    cov_df[col] = cov_df[col] * cov_df["factor"]
    # remove factor
    cov_df = cov_df.drop("factor", axis=1)
    return cov_df


def normalize_GC(cov_df):
    """
    normalize GC for an entire tumor_normal sample
    """
    for col in cov_df.columns:
        if col.startswith("Cov"):
            show_output(f"Normalizing GC ratio for {col}.")
            cov_df = normalize_GC_col(cov_df, col)
    cov_df = cov_df.reset_index(drop=True).sort_values(["Chr", "Pos"])
    return cov_df


def normalize(cov_df):
    """
    normalize general coverage to 100 (no GC normalization)
    """
    for col in cov_df.columns:
        if col.startswith("Cov"):
            show_output(f"Normalizing coverage for {col}.")
            cov_df.loc[:, col] = cov_df[col] / cov_df[col].mean() * 100
    return cov_df


def amazonize(*cov_dfs):
    """
    detect male samples via below-threshold X-chrom coverage
    coverage on chrX is doubled in these samples
    """

    amazon_dfs = []
    for cov_df in cov_dfs:
        # get coverage cols
        cov_cols = [col for col in cov_df.columns if col.startswith("Cov")]
        # create the agg dictionary
        cov_agg = {col: "mean" for col in cov_cols}
        # compute x_coverage for all samples using agg dictionary
        X_coverage = cov_df.query('Chr == "chrX"').agg(cov_agg)
        # filter out the male samples
        male_cols = [col for col in cov_cols if X_coverage[col] < 75]
        # adjust the coverage for male samples
        cov_df.loc[cov_df["Chr"] == "chrX", male_cols] = cov_df[male_cols] * 2
        amazon_dfs.append(cov_df)
    return amazon_dfs


def compute_stats(*dfs):
    """
    get statistics
    """

    compute_dfs = []
    for df in dfs:
        # remove all pre-existing stats
        df = df.drop([col for col in df.columns if col.startswith("PONcov")], axis=1)
        # set index for all non-coverage columns
        index_cols = [col for col in df.columns if not col.startswith("Cov")]
        cov_df = df.drop(index_cols, axis=1)
        df["PONcov_mean"] = cov_df.mean(axis=1)
        df["PONcov_median"] = cov_df.median(axis=1)
        df["PONcov_std"] = cov_df.std(axis=1)
        compute_dfs.append(df)
    return compute_dfs


def remove_outliers(*dfs, std_factor=2.5):
    """
    cycle through all sample cols, remove outliers with difference to PONcov greater than std_factor * std
    """
    filter_dfs = []
    for df in dfs:
        for col in [col for col in df.columns if col.startswith("Cov")]:
            df.loc[
                np.abs(df["PONcov_mean"] - df[col]) / df["PONcov_std"] > std_factor, col
            ] = np.nan
        filter_dfs.append(df)
    return filter_dfs


def make_PON_coverage(
    config={
        "PONcoverage": {
            "stdFactor": 2.5  # only exonPositions straighing within std_factor * std around meanCoverage are kept
        },
        "PON_path": ".",  # path to the PON folder
    },
    save=True,
):
    """ """

    # paths
    pon_path = config["PON_path"]
    # load all sample coverages for one chromosome
    cov_df = combine_PON_coverage(config=config)

    # make chrom categorical
    chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    cov_df.loc[:, "Chr"] = pd.Categorical(cov_df["Chr"], chrom_list)

    # normalize and add mean values and std
    show_output("Normalizing coverage and removing GC dependencies for PON coverage.")
    covGC_df = normalize_GC(cov_df)

    cov_df = normalize(cov_df)

    show_output("Lifting X-coverages for male samples to XX coverage.")
    covGC_df, cov_df = amazonize(covGC_df, cov_df)

    show_output("Computing stats.")
    covGC_df, cov_df = compute_stats(covGC_df, cov_df)

    std_factor = config["PONcoverage"]["stdFactor"]
    show_output("Remove outliers and recompute stats.")
    filterGC_df, filter_df = remove_outliers(covGC_df, cov_df, std_factor=std_factor)

    # save and adjust the output columns
    base_cols = ["Chr", "Pos", "ExonPos"]
    # map_cols = [col for col in cov_df.columns if col.startswith("map")]
    cov_cols = [col for col in cov_df.columns if col.startswith("Cov")]
    stat_cols = [col for col in cov_df.columns if col.startswith("PONcov")]

    cov_df = cov_df.loc[:, base_cols + cov_cols + stat_cols]
    covGC_df = covGC_df.loc[:, base_cols + cov_cols + stat_cols]
    filter_df = filter_df.loc[:, base_cols + stat_cols]
    filterGC_df = filterGC_df.loc[:, base_cols + stat_cols]

    # save dataframes
    if save:
        PON_cov_file = os.path.join(pon_path, "CNV/pon.cov.full.gz")
        show_output(f"Saving combined PON coverage file {PON_cov_file}.")
        cov_df.to_csv(PON_cov_file, sep="\t", index=False, compression="gzip")
        # GC variant
        PON_cov_file = PON_cov_file.replace("full", "fullGC")
        show_output(f"Saving combined PON coverage file {PON_cov_file}.")
        covGC_df.to_csv(PON_cov_file, sep="\t", index=False, compression="gzip")

        # filtered
        PON_cov_file = os.path.join(pon_path, "CNV/pon.cov.filter.gz")
        show_output(f"Saving filtered PON coverage file {PON_cov_file}.")
        filter_df.to_csv(PON_cov_file, sep="\t", index=False, compression="gzip")
        PON_cov_file = PON_cov_file.replace("filter", "filterGC")
        show_output(f"Saving filtered PON coverage file {PON_cov_file}.")
        filterGC_df.to_csv(PON_cov_file, sep="\t", index=False, compression="gzip")

        show_output("Finished", color="success")

    return cov_df, filter_df


# ################# SAMPLE ##############################################
def combine_chrom_cnv(sample, config={}):
    """
    combine the coverages for all chroms and add the GC ratio
    """

    cov_dfs = []
    snp_dfs = []
    chrom_list = [f"chr{c + 1}" for c in range(22)] + ["chrX"]
    for chrom in chrom_list:
        cov_file = os.path.join(config["cov_path"], f"{sample}.{chrom}.cov.gz")
        snp_file = os.path.join(config["snp_path"], f"{sample}.{chrom}.snp")
        show_output(f"Loading coverage file {cov_file}")
        cov_df = pd.read_csv(cov_file, sep="\t", compression="gzip")
        show_output(f"Loading snp file {snp_file}")
        snp_df = pd.read_csv(snp_file, sep="\t")
        cov_dfs.append(cov_df)
        snp_dfs.append(snp_df)
    cov_df = pd.concat(cov_dfs).reset_index(drop=True)
    snp_df = pd.concat(snp_dfs).reset_index(drop=True)
    # make chrom categorical
    cov_df.loc[:, "Chr"] = pd.Categorical(cov_df["Chr"], chrom_list)
    snp_df.loc[:, "Chr"] = pd.Categorical(snp_df["Chr"], chrom_list)
    return cov_df, snp_df


def get_full_exon_pos(*dfs, pon_cov_df):
    """
    adds the accumulated exonic position (over all chroms from PON data) to provided dfs
    """

    # create chrom_df from pon_df
    chrom_df = pon_cov_df.groupby("Chr").agg(dict(ExonPos=["min", "max"]))["ExonPos"]
    chrom_df["chrAdd"] = chrom_df["max"].cumsum().shift(fill_value=0)
    chrom_df = chrom_df.loc[:, "chrAdd"].reset_index()

    full_dfs = []
    for df in dfs:
        # merge with chrom_df
        df = df.merge(chrom_df)
        # get FullExonPos from ExonPos and chrAdd
        df.loc[:, "FullExonPos"] = df["ExonPos"] + df["chrAdd"]

        # get col index of "ExonPos" for adjacent inserting of "FullExonPos"
        insert_index = df.columns.get_loc("ExonPos") + 1
        cols = list(df.columns)
        out_cols = cols[:insert_index] + ["FullExonPos"] + cols[insert_index:-2]
        full_dfs.append(df.loc[:, out_cols])
    return full_dfs


def log2ratio(df, cov_col, pon_cov_col="PONcov_mean"):
    """
    add log2ratio (log2(COV/PONCOV) to coverage data
    """
    # get the incoming columns
    cols = list(df.columns)
    # mask rows where logging does not compute
    loggable = df[cov_col] * df["PONcov_mean"] != 0
    # apply the log
    log_col = cov_col.replace("Cov", "log2ratio")
    df.loc[loggable, log_col] = np.log2(
        df.loc[loggable, cov_col] / df.loc[loggable, pon_cov_col]
    )
    # get col index of cov_col for inserting log_col
    insert_index = df.columns.get_loc(cov_col) + 1
    out_cols = cols[:insert_index] + [log_col] + cols[insert_index:]
    return df.loc[:, out_cols]


def filter_cov(cov_df, config={}):
    """
    takes the config and applies pre-filtering for rolling computation
    """

    c = config["filter"]

    # coverage
    ccov = c["cov"]
    minGC, maxGC = ccov["GCrange"]
    minPONcov = ccov["minPONcov"]
    maxPONstd = ccov["maxPONstd"]
    cov_df = cov_df.query(
        "(@minGC < GCratio < @maxGC) and PONcov_mean >= @minPONcov and PONcov_std < @maxPONstd"
    )
    return cov_df


def filter_snp(snp_df, config={}):
    """
    takes the config and applies pre-filtering for rolling computation
    """

    c = config["filter"]

    csnp = c["snp"]
    # minDepth = csnp['minDepth']
    # minVAF = csnp['minVAF']
    maxPONVAF = csnp["maxPONVAF"]
    minPONDepth = csnp["minPONDepth"]

    snp_df = snp_df.query("@minPONDepth <= PONDepth and PONVAF < @maxPONVAF")

    return snp_df


def combine_sample_CNV(sample, config={}):
    """
    combined snp and cov sample data
    adds FullExonPos to both files

    COV:
        performs GC normalization
        adds the PON coverage and computes log2ratio coverage
    SNP:
        adds PON VAF to sample SNP data for filtering

    """

    # ## LOAD
    show_output(f"Loading data for sample {sample}")
    # combine chrom SNP and COV data for that sample
    cov_df, snp_df = combine_chrom_cnv(sample, config=config)
    show_output("Finished", color="success")

    # ## NORMALIZE
    # get config
    if config["coverage"]["GCnormalize"]:
        show_output("Normalizing coverage using GCratio segmentation")
        cov_df = normalize_GC(cov_df)
        # use GC PON
        PON_cov_file = os.path.join(config["PON_path"], "CNV/pon.cov.filterGC.gz")
    else:
        cov_df = normalize(cov_df)
        # use nonGC PON
        PON_cov_file = os.path.join(config["PON_path"], "CNV/pon.cov.filter.gz")

    # ## INCLUDE PON
    if os.path.isfile(PON_cov_file):
        show_output(f"Loading PON coverage from {PON_cov_file}")
        pon_cov_df = pon_cov_df = pd.read_csv(
            PON_cov_file, sep="\t", compression="gzip"
        )
    else:
        show_output(f"PON coverage file {PON_cov_file} not found!", color="warning")
        return

    show_output("Adding FullExonPos from PON coverage to coverage and snp data.")
    cov_df, snp_df = get_full_exon_pos(cov_df, snp_df, pon_cov_df=pon_cov_df)
    show_output(
        "Merging coverage with PON coverage and performing log2ratio computation."
    )
    cov_df = cov_df.merge(pon_cov_df)

    del pon_cov_df
    for col in cov_df.columns:
        if col.startswith("Cov"):
            cov_df = log2ratio(cov_df, col)
    show_output("log2ratio computation finished.", color="success")

    # PON snp
    PON_snp_file = os.path.join(config["PON_path"], "CNV/pon.snp.gz")
    if os.path.isfile(PON_snp_file):
        show_output(
            f"Loading PON SNP from {PON_cov_file} and merging into sample SNP data."
        )
        pon_snp_df = pon_snp_df = pd.read_csv(
            PON_snp_file, sep="\t", compression="gzip"
        )
    else:
        show_output(f"PON snp file {PON_snp_file} not found!", color="warning")
        return
    snp_df = snp_df.merge(pon_snp_df, how="left")
    show_output(f"Finished combining data for sample {sample}")
    return cov_df, snp_df
