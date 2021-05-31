# try with pd merge
import os
import numpy as np
import pandas as pd
from script_utils_CNV import show_output


def combine_PON_coverage(config={}):
    '''
    combine the PON coverages for all chroms
    '''

    # paths
    pon_path = config['PON_path']

    chrom_list = [f"chr{c + 1}" for c in range(22)] + ['chrX']
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
    cov_df.loc[:, "Chr"] = pd.Categorical(cov_df['Chr'], chrom_list)
    return cov_df


def normalize_GC_col(cov_df, col):
    '''
    normalizes one coverage column for GC ratio
    '''
    # compute the normalizer df
    # for each GCratio, norm_df has the difference of the respective mean from arbitrary norm coverage 100
    # remove chrX for the normalization or male genomes will have slightly greater mean
    norm_df = (100 / cov_df.query('Chr != "chrX"').groupby("GCratio").agg({col: 'mean'})).reset_index().rename({col: "factor"}, axis=1)

    # merge to get the factor
    cov_df = cov_df.merge(norm_df)
    # adjust coverage using the factor from norm_df
    cov_df[col] = cov_df[col] * cov_df['factor']
    # remove factor
    cov_df = cov_df.drop("factor", axis=1)
    return cov_df


def normalize_GC(cov_df):
    '''
    normalize GC for an entire tumor_normal sample
    '''
    for col in cov_df.columns:
        if col.startswith("Cov"):
            show_output(f"Normalizing GC ratio for {col}.")
            cov_df = normalize_GC_col(cov_df, col)
    cov_df = cov_df.reset_index(drop=True).sort_values(['Chr', 'Pos'])
    return cov_df


def amazonize(cov_df):
    '''
    detect male samples via below-threshold X-chrom coverage
    coverage on chrX is doubled in these samples
    '''

    # get coverage cols
    cov_cols = [col for col in cov_df.columns if col.startswith("Cov")]
    # create the agg dictionary
    cov_agg = {col: "mean" for col in cov_cols}
    # compute x_coverage for all samples using agg dictionary
    X_coverage = cov_df.query('Chr == "chrX"').agg(cov_agg)
    # filter out the male samples
    male_cols = [col for col in cov_cols if X_coverage[col] < 75]
    # adjust the coverage for male samples
    cov_df.loc[cov_df['Chr'] == "chrX", male_cols] = cov_df[male_cols] * 2
    return cov_df


def compute_stats(df):
    '''
    get statistics
    '''
    # remove all pre-existing stats
    df = df.drop([col for col in df.columns if col.startswith("PONcov")], axis=1)
    # set index for all non-coverage columns
    index_cols = [col for col in df.columns if not col.startswith("Cov")]
    cov_df = df.drop(index_cols, axis=1)
    df['PONcov_mean'] = cov_df.mean(axis=1)
    df['PONcov_median'] = cov_df.median(axis=1)
    df['PONcov_std'] = cov_df.std(axis=1)
    return df


def remove_outliers(df, std_factor=2.5):
    '''
    cycle through all sample cols, remove outliers with difference to PONcov greater than std_factor * std
    '''
    for col in [col for col in df.columns if col.startswith("Cov")]:
        df.loc[np.abs(df['PONcov_mean'] - df[col]) / df['PONcov_std'] > std_factor, col] = np.nan
    return df


def make_PON_coverage(config={
    'PONcoverage': {
        'stdFactor': 2.5  # only exonPositions straighing within std_factor * std around meanCoverage are kept
    },
    'PON_path': '.',  # path to the PON folder
}):
    '''

    '''
    # load all sample coverages for one chromosome
    cov_df = combine_PON_coverage(config=config)

    # normalize and add mean values and std
    show_output("Normalizing coverage and removing GC dependencies for PON coverage.")
    cov_df = normalize_GC(cov_df)

    show_output("Lifting X-coverages for male samples to XX coverage.")
    cov_df = amazonize(cov_df)

    show_output("Computing stats.")
    cov_df = compute_stats(cov_df)

    std_factor = config['PONcoverage']['stdFactor']
    show_output("Remove outliers and recompute stats.")
    filter_df = remove_outliers(cov_df, std_factor=std_factor)

    # save and adjust the output columns
    base_cols = ['Chr', 'Pos', 'ExonPos']
    # map_cols = [col for col in cov_df.columns if col.startswith("map")]
    cov_cols = [col for col in cov_df.columns if col.startswith("Cov")]
    stat_cols = [col for col in cov_df.columns if col.startswith("PONcov")]

    cov_df = cov_df.loc[:, base_cols + cov_cols + stat_cols]
    filter_df = filter_df.loc[:, base_cols + stat_cols]

    return cov_df, filter_df
