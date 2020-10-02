import os
import pandas as pd
import numpy as np
from script_utils import show_output


def combine_coverage(chrom, sample_file_list):
    cov_df = pd.DataFrame(columns=['ExonPos', 'Pos'])
    for sample_file in sample_file_list:
        if not os.path.isfile(sample_file):
            show_output(f"no file: {sample_file}", color="warning")
            continue
        sample = os.path.basename(sample_file).split('.')[0]
        show_output(f"Reading {sample} from {sample_file}.")
        df = pd.read_csv(sample_file, sep='\t', compression='gzip').loc[:, [
            'Pos', 'ExonPos', 'Coverage']].rename(columns={'Coverage': sample})
        cov_df = cov_df.merge(df, on=['ExonPos', 'Pos'], how='outer')
    cov_df = cov_df.fillna(0).sort_values('ExonPos')
    cov_df['Chr'] = chrom
    # reorder columns
    cols = ['Chr', 'Pos', 'ExonPos'] + list(cov_df.columns)[2:-1]
    return cov_df.loc[:, cols]


def normalize_coverage(cov_df, normCov=100):
    norm_df = cov_df.set_index(['Chr', 'Pos', 'ExonPos'])
    norm_df = norm_df / norm_df.mean() * normCov
    return norm_df.reset_index()


def add_mean(norm_df):
    norm_df = norm_df.set_index(['Chr', 'Pos', 'ExonPos'])
    norm_df['meanCov'] = norm_df.mean(axis=1)
    norm_df['medianCov'] = norm_df.median(axis=1)
    norm_df['std'] = norm_df.std(axis=1)
    return norm_df.reset_index()


def filter_coverage(df, minCov=20, maxMeanSTD=20):
    filter_df = df.query('meanCov > @minCov and std < @maxMeanSTD')
    return filter_df


def remove_outliers(df, std_factor=2.5):
    '''
    cycle through all sample cols, remove outliers and recompute the means and std
    '''

    for col in list(df.columns)[3:-3]:
        df.loc[np.abs(df['meanCov'] - df[col]) / df['std']
               > std_factor, col] = np.nan
    return add_mean(df.iloc[:, :-3])


def make_PON_coverage(chrom, sample_list, config={
    'normCov': 100,       # to what value are coverages normalized
    'minCov': 20,        # only exonPositions with the average coverage above minCov are kept
    'maxMeanSTD': 20,  # only exonPositions with a coverage std below max_mean_std are kept
    # only exonPositions straighing within stdFactor * std around meanCoverage are kept
    'stdFactor': 3,
}):
    '''
    load the coverages for all the PON files for one chromosome and write that to a file if given
    '''
    # load all sample coverages for one chromosome
    cov_df = combine_coverage(chrom, sample_list)

    # normalize and add mean values and std
    mean_df = add_mean(normalize_coverage(cov_df, normCov=config['normCov']))

    # filter hard regions and outlying data points
    filter_df = filter_coverage(
        mean_df, minCov=config['minCov'], maxMeanSTD=config['maxMeanSTD'])
    filter_removed_df = remove_outliers(
        filter_df, std_factor=config['stdFactor'])

    # output only if outpath is given

    return mean_df, filter_removed_df
