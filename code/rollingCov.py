import pandas as pd

from rollingCNV import llh, rolling_data, interpolate
from script_utils_CNV import show_output


def compute_covLLH_col(df, log_col, params={}):
    """
    computes the local log-likelihood of belonging to the center gaussian for specified log_column
    """

    cols = list(df.columns)

    min_log2ratio, max_log2ratio = params["center_range"]
    # get the sigma and mean of the center band log2ratio
    center_logs = df.query(f"@min_log2ratio < {log_col} < @max_log2ratio")[log_col]
    sigma = center_logs.std() * params["sigma_factor"]
    mean = center_logs.mean()
    show_output(
        f"Computing log-likelihood of {log_col} belonging to center gaussian [mean:{round(mean, 3)}, sigma:{round(sigma,3)}]"
    )
    llh_col = log_col.replace("log2ratio", "covLLH")
    df.loc[:, llh_col] = llh(df[log_col], mean=mean, sigma=sigma)
    # get col index of log_col for adjacent inserting of llh_col

    insert_index = df.columns.get_loc(log_col) + 1
    df = df.loc[:, cols[:insert_index] + [llh_col] + cols[insert_index:]]
    return df


def compute_covLLH(df, config={}):
    """
    computes the local log-likelihood of belonging to the center gaussian
    """

    for col in df.columns:
        if col.startswith("log2ratio"):
            df = compute_covLLH_col(df, col, params=config["rolling"]["cov"]["LLH"])

    return df


def rolling_coverage(cov_df, config={}):
    """
    cycle through the chroms and perform rolling window computations of data set in config
    """

    # split the params dict for easier access
    cc = config["rolling"]["cov"]
    data_params = cc["data"]

    # perform llh computation for coverage if needed
    if "covLLH" in data_params:
        cov_df = compute_covLLH(cov_df, config=config)

    cov_df = rolling_data(cov_df, data_params=data_params, roll_config=cc)

    return cov_df


def interpolate_fullexonpon(merge_df):
    chrom_dfs = []
    for chrom in merge_df["Chr"].unique():
        chrom_df = merge_df.query("Chr == @chrom")
        chrom_df = interpolate(
            chrom_df, "FullExonPos", ref_col="Pos", expand_limit=1000000
        )
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values("FullExonPos")
    df.loc[:, "FullExonPos"] = df["FullExonPos"].fillna(0).astype(int)
    return df
