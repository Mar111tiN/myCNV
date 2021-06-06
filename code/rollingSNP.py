import pandas as pd
import numpy as np
import re

from rollingCNV import rolling_data, one_col_rolling, interpolate
from script_utils_CNV import show_output


def mergeCov2SNP(snp_df, cov_df):
    """
    merge the columns to get all the coverage data from the rolling_coverage
    """

    # get the columns
    snp_cols = list(snp_df.columns)
    base_cols = snp_cols[:4]
    cov_cols = [
        col
        for col in cov_df.columns
        if col.endswith("mean") or col.endswith("sum") or col == "GCratio"
    ]
    # reduce cols of cov_df
    cov_df = cov_df.loc[:, base_cols + cov_cols]
    # PONVAF and PONDepth have to be filled to zero
    for col in ["PONVAF", "PONDepth"]:
        snp_df.loc[:, col] = snp_df[col].fillna(0)

    snp_df = (
        snp_df.merge(cov_df, how="outer")
        .sort_values("FullExonPos")
        .reset_index(drop=True)
    )

    # interpolate the data
    for col in cov_cols:
        show_output(f"Interpolating {col}")
        snp_df = interpolate(snp_df, col, expand_limit=200)

    # reduce to VAF values
    snp_df = snp_df.query("PONVAF == PONVAF")
    return snp_df


def center_vaf(snp_df):
    """
    adjusts VAFs to 0.5
    also includes the off-center VAF offVAF and
    absVAF
    """

    df = snp_df.copy()
    cols = list(df.columns)

    for col in df.columns:
        if col.startswith("VAF"):
            cols = list(df.columns)
            insert_index = df.columns.get_loc(col) + 1
            df.loc[:, col] = df[col] + 0.5 - df[col].mean()
            # get additional features from VAFs
            # get offVAF as measure of straighing from center
            df.loc[:, f"off{col}"] = (df[col] - 0.5) * 2
            # get absVAF as absolute measure of straighing from center
            df.loc[:, f"abs{col}"] = np.abs(df[f"off{col}"])
            out_cols = (
                cols[:insert_index] + [f"abs{col}", f"off{col}"] + cols[insert_index:]
            )
            df = df.loc[:, out_cols]
    return df


# get the density computer for rolling
def make_get_density(window_size=20):
    """
    helper for returning a density computer for given window_size
    """

    def SNPdensity(data):
        return window_size / (data.max() - data.min())

    return SNPdensity


def get_fallSNP(snp_df, config={}):
    """
    computes snpDensity and cumulative metrix snpFall (SNPdensity * PONVAF)
    """
    df = snp_df.copy()

    c = config["rolling"]["snp"]

    window = c["data"]["SNPdensity"]
    # create the callback for rolling window computation
    get_SNPdensity = make_get_density(window)

    # make params for rolling_data function and roll custom function
    fallSNP_params = {"FullExonPos": {get_SNPdensity: window}}
    df = rolling_data(df, data_params=fallSNP_params, roll_config=c)
    # compute fallSNP from density
    df.loc[:, "fallSNP"] = df["SNPdensity"] * df["PONVAF"]
    # normalize
    df.loc[:, "fallSNP"] = df["fallSNP"] / df["fallSNP"].max()
    return df


def filter_snp(snp_df, config={}):
    """
    takes the config and applies pre-filtering for rolling computation
    """

    df = snp_df.copy()
    c = config["filter"]["snp"]

    # ### build the query

    # minDepth = csnp['minDepth']
    # minVAF = csnp['minVAF']

    # PONVAF and PONDepth have to be filled to zero
    for col in ["PONVAF", "PONDepth"]:
        df.loc[:, col] = df[col].fillna(0)
    # pon query
    # minoffVAF = c["minoffVAF"]
    PON_query = f"PONVAF < {c['maxPONVAF']}"
    # Fall_query = f"(Fall2 < {maxFall} and offVAF2_sum > {minoffVAF})"
    Fall_query = f"fallSNP < {c['maxFallSNP']} and SNPdensity < {c['maxSNPdensity']}"
    # map query
    # extract map query from filter config
    # 'map30_0 > 0.1 and map50_0 > 0.1 and map75_1 > 0.1 and map100_2 > 0.1'
    map_query = " and ".join(
        [f"{m} >= {c[m]}" for m in c.keys() if m.startswith("map")]
    )

    snp_query = f"{PON_query} and {map_query} and {Fall_query}"
    show_output(f'Filtering SNP data using {snp_query.replace("@", "")}')
    df = df.query(snp_query)

    return df


def rolling_snp(snp_df, roll_cov_df, config={}):
    """
    performs data preprocessing for CNV analysis
    """

    df = snp_df.copy()
    # get config

    show_output("Starting rolling SNP processing.", time=True)

    show_output("Merging SNP and coverage data.")
    df = mergeCov2SNP(df, roll_cov_df)

    show_output("Centering VAFs.")
    df = center_vaf(df)

    show_output("Running fallSNP detection.")
    # compute fallSNP
    show_output("Removing fall SNPs.")
    df = get_fallSNP(df, config=config)
    # filter snp

    show_output("Reducing noise...")
    df = filter_snp(df, config=config)
    show_output("Finished SNP processing!", time=True, color="success")
    return df


def remergeCNV(snp_df, cov_df):
    """
    remerge the rolling data for coverage and snp for an integrated data file
    containing many nas
    """

    # select the right cols

    base_cols = ["Chr", "Pos", "ExonPos", "FullExonPos"]
    snp_cols = [col for col in snp_df.columns if col.startswith("VAF")]
    log2_pat = re.compile(r"log2ratio[0-9]+(_mean)?$")
    cov_cols = [col for col in cov_df.columns if re.match(log2_pat, col)]
    snp_df = snp_df.loc[:, base_cols + snp_cols]
    cov_df = cov_df.loc[:, base_cols + cov_cols]
    cnv_df = (
        snp_df.merge(cov_df, how="outer")
        .sort_values("FullExonPos")
        .reset_index(drop=True)
    )
    return cnv_df


# ############# LEGACY ###############


# def compute_snp_llh(df, mean=0.5, sigma=0.2):
#     """
#     computes the local log-likelihood of belonging to the center gaussian
#     """

#     show_output(
#         f"Computing log-likelihood of VAF belonging to center gaussian [mean:{round(mean, 3)}, sigma:{round(sigma,3)}]"
#     )
#     df.loc[:, "snpLLH"] = llh(df["VAF"], mean, sigma)

#     # for homoSNPs reduce the VAFs to the ones above mean
#     upper_vafs = df.query("@mean < VAF")["VAF"]
#     # then compute the hsnpLLH

#     show_output(
#         f"Computing log-likelihood of VAF belonging to purity100  [mean:1, sigma:{round(sigma,3)}]"
#     )
#     # these are called hsnp
#     # upper_vafs only contains half the snps, the remaining have to be interpolated
#     df.loc[:, "hsnpLLH"] = llh(upper_vafs, 1, sigma)
#     df = interpolate(df, "hsnpLLH", expand_limit=50)
#     return df


# def remove_fallSNP(snp_df, mean=0.5, std=0.2, params={}):
#     """
#     removes the falling SNP probably caused by mismapping
#     """

#     window = params["offVAFwindow"]
#     cutoff = params["maxFallSNP"]

#     # get the density computer for rolling
#     get_SNPdensity = make_get_density(window)
#     # cycle through chroms
#     chrom_dfs = []
#     for chrom in snp_df["Chr"].unique():
#         df = snp_df.query("Chr == @chrom").reset_index(drop=True)
#         if len(df.index) < 10:
#             continue
#         # get the snp
#         df = one_col_rolling(
#             df,
#             df.query("VAF < 0.95"),
#             "ExonPos",
#             get_SNPdensity,
#             window_size=window,
#             diff_exp=4,
#         )
#         df.loc[:, "SNPdensity"] = df["SNPdensity"] / df["SNPdensity"].mean()

#         # get the offVAFsum
#         df = one_col_rolling(
#             df,
#             df.query("VAF < 0.95"),
#             "offVAF",
#             "sum",
#             window_size=window,
#             normalize=True,
#             diff_exp=4,
#         )

#         # combine both metrices
#         df.loc[:, "fallSNP"] = df["SNPdensity"] * df["offVAFsum"]
#         # now remove the ones below average VAFstd
#         df = df.query("VAF > @mean - @std / 2 or fallSNP > @cutoff")
#         chrom_dfs.append(df)

#     return pd.concat(chrom_dfs).sort_values("FullExonPos").reset_index(drop=True)


# def expand_SNPdata(snp_df, config):
#     """
#     retrieve a few data columns locally to use rolling windows on
#     this needs to be done chromosome-wise in order to avoid gap effects
#     VAF limits are also applied here
#     """

#     # split the params dict for easier access
#     params = config["snp"]
#     filter_params = params["filter"]
#     # data_params = params['data']

#     # reduce the snp_df using lower config limit
#     # upper limit has to be set later as we still need the homoSNP llh
#     VAFmin, VAFmax = filter_params["VAF"]
#     snp_df = snp_df.query("@VAFmin < VAF")

#     # get std and mean of VAF
#     minVAF, maxVAF = params["LLH"]["center_range"]
#     # get the sigma and mean of the center band VAF (extracted as pd.Series center_vafs)
#     center_vafs = snp_df.query("@minVAF < VAF < @maxVAF")["VAF"]
#     # get width of gaussian from std * sigma_factor
#     VAFstd = center_vafs.std()
#     VAFmean = center_vafs.mean()

#     # get additional features from VAFs
#     snp_df.loc[:, "offVAF"] = (snp_df["VAF"] - VAFmean) * 2
#     # absolute values for cluster
#     snp_df.loc[:, "absVAF"] = np.abs(snp_df["offVAF"])

#     ########## remove fallSNP ########
#     fs_params = params["fallSNP"]
#     if fs_params["run"]:
#         show_output("Removing falling SNPs")
#         snp_df = remove_fallSNP(snp_df, mean=VAFmean, std=VAFstd, params=fs_params)

#     ######## LLH  #####################
#     # get the snpLLH and hsnpLLH
#     # get config params
#     sigma = VAFstd * params["LLH"]["sigma_factor"]
#     # hsnpLLH is computed in order to rescue high absVAF that would have been filtered out
#     # lower VAF is already removed because density of VAF ~0 is highly irregular and would confound density estimates
#     snp_df = compute_snp_llh(snp_df, mean=VAFmean, sigma=sigma)
#     return snp_df.query("VAF < @VAFmax").reset_index(drop=True)


# def rolling_SNP(snp_df, config):
#     """
#     cycle through the chroms and perform rolling window computations of snp data set in config
#     """

#     # split the params dict for easier access
#     params = config["snp"]
#     filter_params = params["filter"]
#     data_params = params["rolling_data"]
#     debug = config["debug"]

#     minDepth = filter_params["minDepth"]
#     filter_df = snp_df.query("Depth >= @minDepth")
#     show_output("Performing rollingSNP computations.")
#     rolling_df = rolling_data(
#         snp_df,
#         filter_df,
#         expand=params["expand"],
#         ddof=config["ddof"],
#         debug=debug,
#         data_params=data_params,
#     )
#     return rolling_df


# def apply_rolling_SNP(snp_df, config):

#     # get extra data
#     snp_df = expand_SNPdata(snp_df, config)
#     # do the rolling
#     snp_df = rolling_SNP(snp_df, config)
#     # get the CNV and Center blocks
#     snp_df = get_CNV_blocks(snp_df, "snpLLH", config)

#     # select columns for output
#     base_cols = list(snp_df.columns[:4])

#     snp_cols = [
#         col
#         for col in snp_df.columns[4:]
#         if not "log2" in col and not "cov" in col and not "off" in col
#     ]
#     rolling_snp_df = snp_df[base_cols + snp_cols]
#     cluster_cols = [
#         "log2ratio",
#         "log2ratiomean",
#         "VAF",
#         "absVAF",
#         "absVAFmean",
#         "snpLLHsum",
#     ]
#     cluster_cols += [col for col in snp_df.columns if "Center" in col or "CNV" in col]
#     cluster_df = snp_df[base_cols + cluster_cols]
#     return rolling_snp_df, cluster_df
