import pandas as pd
import numpy as np
from script_utils_CNV import show_output


def interpolate(df, data_col, ref_col="FullExonPos", expand_limit=20):
    """
    interpolates missing values in data_col using linear interpolation based on ref_col
    """

    cols = list(df.columns)
    # set FullExonPos as index for the interpolation method to work on proper intervals
    df = df.reset_index(drop=False).set_index(ref_col, drop=False)
    df.loc[:, data_col] = df[data_col].interpolate(
        method="values", limit=expand_limit, limit_direction="both"
    )
    return df.set_index("index")[cols]


def normalize_df(df, col):
    """
    normalize a column of a df
    """
    _min = df[col].min()
    _max = df[col].max()
    df.loc[:, col] = (df[col] - _min) / (_max - _min)
    return df


def llh(data, mean, sigma):
    """
    compute the density function for a given gaussian
    takes a pd.Series or np.array
    """
    s = np.sqrt(2 * np.pi) * sigma
    return np.exp((data - mean) ** 2 / (-2 * (sigma ** 2))) / s


def one_col_rolling(
    df,
    df_filter,
    col,
    aggr,
    window_size=200,
    expand_limit=20,
    normalize=False,
    debug=False,
    diff_exp=2,
    ddof=0,
):
    """
    performs rolling computation of <agg> on data column <col> with given window size
    the aggregation can be a:
        - callable taking df[col] as argument and returning a scalar
            column name will be taken from function name (stripping underscores)
        - string expression understood by the agg-function of the pandas.groupby API
            column name will be composed of col + aggr
    computation is performed on a left and right rolling window
    missing margins are filled by the counterpart window function
    a diff column is included ()
    """

    org_cols = list(df.columns)
    # rolling left
    # check if aggr is a function
    if callable(aggr):
        if debug:
            show_output(f"Aggregating custom function {aggr.__name__}")
        df.loc[:, "L"] = df_filter[col].rolling(window_size).apply(aggr)
        # pass the function name for ensuing column naming
        col_name = aggr.__name__.replace("_", "")
    else:
        # get the right computation by passing aggr to .agg()
        # only this allows passing methods as string
        df.loc[:, "L"] = df_filter[col].rolling(window_size).agg(aggr, ddof=ddof)
        col_name = col + aggr

    # rolling right by shifting the L column
    df.loc[:, "R"] = df.shift(-window_size + 1)["L"]

    diff_name = col_name + "Diff"
    new_cols = org_cols + [col_name, diff_name]
    if debug:
        new_cols += [f"{col_name}L", f"{col_name}R"]
    # skips interpolation if value == 0
    if interpolate:
        # interpolate missing values
        for c in ["L", "R"]:
            df = interpolate(df, c, expand_limit=expand_limit)
    # fill the margins
    L_margin = df["L"].first_valid_index()
    df.loc[:L_margin, "L"] = df["R"]

    ######################
    # this fails if window is too large!!!!
    # debug
    #####################

    R_margin = df["R"].last_valid_index() + 1
    try:
        df.loc[R_margin:, "R"] = df["L"]
    except:
        R_margin -= 1
        ##### DEBUG ###########
        print("<DEBUG> ADJUSTING R_MARGIN")
        #######################
        df.loc[R_margin:, "R"] = df["L"]

    # get the Diff
    df.loc[:, diff_name] = np.abs(df["R"] - df["L"])
    # normalize to max
    df.loc[:, diff_name] = df[diff_name] / df[diff_name].max()
    # here, contribution of L and R is controlled by diff value
    df.loc[:, col_name] = df["R"] * df[diff_name] + df["L"] * (1 - df[diff_name])

    if normalize:
        df = normalize_df(df, col_name)
        if debug:
            for c in ["L", "R"]:
                df = normalize(df, c)

    # square the diff
    df.loc[:, diff_name] = df[diff_name] ** diff_exp

    if debug:
        # specify col names of L and R
        df = df.rename(columns=dict(L=f"{col_name}L", R=f"{col_name}R"))

    # reduce to the right columns
    return df.loc[:, new_cols]


def rolling_data(df, filter_df, expand=0.25, ddof=0, debug=False, data_params={}):
    """
    cycles through the data params (rolling_data object from config dict)
    and performs rolling computations for these params
    generic function to be used by rolling_coverage and rolling_SNP
    """

    # now do global normalization for sum aggregations:
    # cycle through rolling_data
    for data_col in data_params.keys():
        for agg in data_params[data_col].keys():
            # cycle through the chroms
            chrom_dfs = []
            for chrom in df["Chr"].unique():
                # get the chrom_dfs
                chrom_df = df.query("Chr == @chrom").sort_values("FullExonPos")
                filter_chrom_df = filter_df.query("Chr == @chrom").sort_values(
                    "FullExonPos"
                )
                window_size = data_params[data_col][agg]
                if len(filter_chrom_df.index) < 2 * window_size:
                    continue
                expand_limit = int(expand * window_size)
                # show_output(
                #     f"Computing rolling window for {agg} of {data_col} with window size {window_size} on {chrom}"
                # )
                chrom_df = one_col_rolling(
                    chrom_df,
                    filter_chrom_df,
                    data_col,
                    agg,
                    window_size=window_size,
                    expand_limit=expand_limit,
                    ddof=ddof,
                )
                chrom_dfs.append(chrom_df)
            # combine the chrom_dfs
            df = pd.concat(chrom_dfs).sort_values("FullExonPos")

            # Normalization
            # only do normalization for sum aggregations
            if not agg == "sum":
                continue
            print(f"Normalizing {data_col} {agg}")
            # get the columns for normalization
            col_name = data_col + agg
            cols = [col_name]
            if debug:
                cols += [f"{col_name}L", f"{col_name}R"]
            for c in cols:
                _min = df[c].min()
                _max = df[c].max()
                df.loc[:, c] = (df[c] - _min) / (_max - _min)
    return df


def get_blocks(df, col, min_size=0):
    """
    takes a column of binary containment to certain group and returns block numbers and respective block_sizes
    excludes blocks below a certain block size limit
    """

    org_cols = list(df.columns)
    # find gaps where col value changes
    df.loc[:, ["gap"]] = (df[col] != df.shift(1)[col]).astype(int)
    # set the col names
    blocksize = f"{col}block_size"
    block = f"{col}block"

    # enumerate the gaps for blockID where value is 1
    df.loc[:, [block]] = df["gap"].cumsum() * df[col]
    # group by blocks and count size
    blocks = df.groupby(block)["gap"].count().rename(blocksize)
    # merge block size into df
    df = df.merge(blocks, left_on=block, right_index=True)
    # remove miniscule blocks
    df.loc[df[blocksize] < min_size, block] = 0

    # maybe adjust block labelling to be numerically ordered
    # should maybe done
    df.loc[:, col] = df[block]
    # cols = org_cols + [block, blocksize]
    return df[org_cols]


def get_CNV_blocks(df, data, config):
    """
    finds blocks of CNV with center LLH below threshold
    then it reduces these blocks to the regions within the Diff-peaks..
    ..to only enter the most meaningful data into clustering

    do the same thing for covCenter with values above center threshold
    (definitely belonging to the center)
    """

    # extract params from config
    # for covLLH --> lookup combine.cov.LLH_cutoff
    t = data.replace("LLH", "")
    params = config[t]

    LLH_params = params["LLH_cutoff"]

    col = data + "sum"

    # get boolint whether LLH falls below threshold
    # fillna(1) to exclude any missing coverages
    df.loc[:, f"{t}CNV"] = (df[col].fillna(1) < LLH_params["cnv"]).astype(int)
    # get the covCNV blocks
    df = get_blocks(df, f"{t}CNV", min_size=LLH_params["min_block_size"])
    # reduce data to within Diff-peaks
    df.loc[:, f"{t}CNVcore"] = (
        (df[f"{t}CNV"] > 0) & (df[f"{t}LLHsumDiff"] < LLH_params["max_diff"])
    ).astype(int)

    # here I could also expand the covCNV to the peak of the Diff for covCNVexp
    # get the window_size for core expanding
    window_size = params["rolling_data"][data]["sum"]

    # get the covCenter
    df.loc[:, f"{t}Center"] = (df[col].fillna(0) > LLH_params["center"]).astype(int)
    # get the covCenter blocks
    df = get_blocks(df, f"{t}Center", min_size=LLH_params["min_block_size"])
    # reduce data to within Diff-peaks
    df.loc[:, f"{t}Centercore"] = (
        (df[f"{t}Center"] > 0) & (df[f"{t}LLHsumDiff"] < LLH_params["max_diff"])
    ).astype(int)

    # core expanding
    # get the window_size for core expanding
    window_size = params["rolling_data"][data]["sum"]
    diff_col = data + "Diff"
    # get the covCNV

    return df


####################### LEGACY #############################################
def rollingCNV(sample, sample_cnv_path, pon_cov_path, config):
    """
    combines all the hetSNP and coverage data per sample and
    performs rolling computations for clustering
    returns the combined raw data and the (optionally na_removed) rolling data
    """

    # combine the chromosome data and associate coverage data with pon coverage
    snp_df, cov_df = get_covNsnp(
        sample,
        sample_cnv_path=sample_cnv_path,
        pon_cov_path=pon_cov_path,
        verbose=config["debug"],
    )

    # apply rolling coverage
    show_output(f"Performing rolling coverage computation for sample {sample}.")
    snpcov_df, rolling_cov_df = apply_rolling_coverage(snp_df, cov_df, config)

    # apply rolling SNP
    show_output(f"Performing rolling computation for hetSNP data of sample {sample}.")
    rolling_snpcov_df, cluster_df = apply_rolling_SNP(snpcov_df, config)
    show_output(f"Finished computations for sample {sample}.")
    if config["na_remove"]:
        rolling_cov_df = rolling_cov_df.dropna()
        rolling_snpcov_df = rolling_snpcov_df.dropna()
        show_output(f"Removed missing values from rolling data for sample {sample}.")
    return cov_df, snp_df, rolling_cov_df, rolling_snpcov_df, cluster_df
