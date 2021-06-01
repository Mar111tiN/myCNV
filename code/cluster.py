import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture as GMM


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


def get_centers(merge_df, runs=25, comps=3, VAF_limits=(0.05, 0.95), exclude_X=True):
    """
    use GMM to identify the center cluster and get the means from that
    because GMM occasionally does not identify the center cluster,
    I let the GMM proceed several times and minimize the center cluster
    next, the center cluster can be identified as the maximum center
    """
    VAFmin, VAFmax = VAF_limits
    # fit the centers to the data
    if exclude_X:
        merge_df = merge_df.query('Chr != "chrX"')
    X = merge_df.query("@VAFmin < VAF < @VAFmax and log2ratiomean == log2ratiomean")[
        ["log2ratiomean", "VAF"]
    ]

    gmm = GMM(n_components=comps, covariance_type="diag", n_init=runs).fit(X)
    labels = gmm.predict(X)
    # get the size of the
    _, counts = np.unique(labels, return_counts=True)
    maxcount = np.max(counts)
    centers = pd.DataFrame(gmm.means_, columns=["log2ratio", "VAF"])
    # get mean_cov and meanVAF from largest cluster
    meanCov, meanVAF = centers.loc[np.argmax(counts)]
    size = maxcount

    print(
        f"GMM using {runs} inits: center size {size} meanVAF = {round(meanVAF, 2)} meanCov={round(meanCov, 2)}"
    )

    return meanCov, meanVAF, centers


def center_data(snp_df, config):
    """
    retrieve the centers for scaling using GMM
    """

    meanCov, meanVAF, _ = get_centers(
        snp_df, VAF_limits=config["heteroSNP"]["filter"]["VAF"]
    )
    # center coverage
    if config["coverage"]["center"]:
        print("log2ratio centered around", meanCov)
        snp_df.loc[:, "log2ratiomean"] = snp_df["log2ratiomean"] - meanCov
    if config["heteroSNP"]["center"]:
        print("heteroSNP centered around", meanVAF)
        snp_df.loc[:, "VAF"] = snp_df["VAF"] - meanVAF + 0.5
    return snp_df
