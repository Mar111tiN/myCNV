import os
import pandas as pd
import numpy as np
from script_utils import show_output


def gather_PONcoverage_chrom(sample_df, config={}):
    cov_df = pd.DataFrame(columns=["ExonPos", "Pos"])
    for _, row in sample_df.iterrows():
        sample_file = row["file"]
        sample_name = row["sample"]
        if not os.path.isfile(sample_file):
            show_output(f"no file: {sample_file}", color="warning")
            continue
        if config["verbose_output"]:
            show_output(f"Reading {sample_name} from {sample_file}.")
        df = (
            pd.read_csv(sample_file, sep="\t", compression="gzip")
            .loc[:, ["Pos", "ExonPos", "Coverage"]]
            .rename(columns={"Coverage": sample_name})
        )
        cov_df = cov_df.merge(df, on=["ExonPos", "Pos"], how="outer")
    cov_df = cov_df.fillna(0).sort_values("ExonPos")
    cov_df["Chr"] = sample_df["Chr"].iloc[0]
    # reorder columns
    cols = ["Chr", "Pos", "ExonPos"] + list(cov_df.columns)[2:-1]
    return cov_df.loc[:, cols]


def gather_PONcoverage(sample_df, config={}):
    """
    combine the PONcoverage for all chromosomes
    """

    cov_dfs = []
    for chrom in sample_df["Chr"].unique():
        show_output(f"Collecting PON coverages for {chrom}")
        cov_df = gather_PONcoverage_chrom(
            sample_df.query("Chr == @chrom"), config=config
        )
        cov_dfs.append(cov_df)
    cov_df_full = pd.concat(cov_dfs).reset_index(drop=True)
    return cov_df_full


def normalize_coverage(cov_df, norm_cov=100):
    norm_df = cov_df.set_index(["Chr", "Pos", "ExonPos"])
    norm_df = norm_df / norm_df.mean() * norm_cov
    return norm_df.reset_index()


def equalize_X(norm_df):
    """
    detects samples with a normed chromX coverage below 75 (XY)
    and doubles respective coverages for an overall diploidy for
    normalization of coverage
    """

    # extract x chrom from normalized df and set index to hide non-coverage cols
    x_df = norm_df.query('Chr == "chrX"').set_index(["Chr", "Pos", "ExonPos"])
    no_x_df = norm_df.query('Chr != "chrX"')
    # double the values for samples with mean below 75
    x_df.loc[:, x_df.mean() < 75] = x_df * 2
    # concat norm_df without x and harmonized X chrom df
    equalX_df = pd.concat([no_x_df, x_df.reset_index()]).sort_values(["Chr", "Pos"])
    return equalX_df


def add_mean(norm_df):
    norm_df = norm_df.set_index(["Chr", "Pos", "ExonPos"])
    norm_df["meanCov"] = norm_df.mean(axis=1)
    norm_df["medianCov"] = norm_df.median(axis=1)
    norm_df["std"] = norm_df.std(axis=1)
    return norm_df.reset_index()


def get_full_exon_pos(df):
    """
    adds the accumulated exonic position (over all chroms)
    """

    # save the output columns
    cols = list(df.columns)
    df = df.reset_index(drop=True)
    # adds the last ExonPos of chrom to start of next chromosome
    df.loc[:, "chromStep"] = df.shift(1)["ExonPos"].fillna(0).astype(int)
    df.loc[df["Chr"] == df.shift(1)["Chr"], "chromStep"] = 0
    df["chromAccum"] = df["chromStep"].cumsum()
    df["FullExonPos"] = df["ExonPos"] + df["chromAccum"]
    cols = cols[:2] + ["FullExonPos"] + cols[2:]
    return df[cols]


def remove_outliers(df, std_factor=2.5):
    """
    cycle through all sample cols, remove outliers
    """

    for col in list(df.columns)[3:-3]:
        df.loc[np.abs(df["meanCov"] - df[col]) / df["std"] > std_factor, col] = np.nan
    return add_mean(df.iloc[:, :-3])


def make_PON_coverage(
    sample_df,
    config={
        "normCov": 100,
        # only exonPositions straighing within std_factor * std around meanCoverage are kept
        "stdFactor": 3,
        "verbose_output": False,
    },
):
    # load all sample coverages for one chromosome
    show_output(
        f"Loading all PON coverages from for exome-wide normalization",
        time=True,
    )
    cov_df = gather_PONcoverage(sample_df, config=config)

    # normalize and add mean values and std
    show_output(f"Normalizing PON coverages and adjusting X chrom coverage", time=True)
    eqX_df = equalize_X(normalize_coverage(cov_df, norm_cov=config["normCov"]))
    show_output(f"Re-Normalizing X-adjusted PON coverages", time=True)
    normX_df = normalize_coverage(eqX_df, norm_cov=config["normCov"])
    mean_df = add_mean(normX_df)
    # add full exon coords to normalized PON coverage
    full_df = get_full_exon_pos(mean_df)
    show_output(f"Removing outliers", time=True)
    filter_df = remove_outliers(mean_df, std_factor=config["stdFactor"])
    # remove sample columns and addd full exon coords to filtered PON coverage
    filter_df = get_full_exon_pos(
        filter_df.loc[:, ["Chr", "Pos", "ExonPos", "meanCov", "medianCov", "std"]]
    )
    show_output(f"Finished filtering of PON coverages", time=True, color="success")
    return full_df, filter_df
