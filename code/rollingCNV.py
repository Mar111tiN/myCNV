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


def one_col_rolling(df, col, aggr="mean", window_size=200, roll_config={}):
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
    # UNPACK PARAMS
    normalize = roll_config["normalize"]
    debug = roll_config["debug"]
    diff_exp = roll_config.get("diffexp", 2)
    ddof = roll_config.get("ddof", 0)

    # save org cols
    cols = list(df.columns)

    # ###### ROLLING LEFT
    # check if aggr is a function
    if callable(aggr):
        if debug:
            show_output(f"Aggregating custom function {aggr.__name__}")
        df.loc[:, "L"] = df[col].rolling(window_size).apply(aggr)
        # pass the function name for ensuing column naming
        col_name = aggr.__name__.replace("_", "")
    else:
        # get the right computation by passing aggr to .agg()
        # only this allows passing methods as string
        df.loc[:, "L"] = df[col].rolling(window_size).agg(aggr, ddof=ddof)
        col_name = col + aggr

    # ###### ROLLING RIGHT
    # rolling right by shifting the L column
    df.loc[:, "R"] = df.shift(-window_size + 1)["L"]

    # ###### DIFFING
    diff_name = col_name + "Diff"

    added_cols = [col_name, diff_name]
    if debug:
        added_cols += [f"{col_name}L", f"{col_name}R"]
    # skips interpolation if value == 0
    if interpolate:
        # interpolate missing values
        for c in ["L", "R"]:
            df = interpolate(df, c, expand_limit=10)
    # fill the margins
    try:
        L_margin = df["L"].first_valid_index()
        df.loc[:L_margin, "L"] = df["R"]
        R_margin = df["R"].last_valid_index() + 1
        df.loc[R_margin:, "R"] = df["L"]
    except Exception as e:
        show_output(
            f"An error occurred attempting to fill the margins! {e}", color="warning"
        )

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

    # steepen the diff
    df.loc[:, diff_name] = df[diff_name] ** diff_exp

    if debug:
        # specify col names of L and R
        df = df.rename(columns=dict(L=f"{col_name}L", R=f"{col_name}R"))

    # insert the cols
    # get col index of cov_col for inserting log_col
    insert_index = df.columns.get_loc(col) + 1
    out_cols = cols[:insert_index] + added_cols + cols[insert_index:]

    # reduce to the right columns
    return df.loc[:, out_cols]


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
