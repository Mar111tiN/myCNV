import os
import pandas as pd


def combine_coverage(PONcov_path, chrom, sample_list):
    cov_df = pd.DataFrame(columns=['ExonPos', 'Pos'])
    for sample in sample_list:
        file = os.path.join(PONcov_path, f"{sample}.{chrom}.bedCov")
        if not os.path.isfile(file):
            continue
        print(f"Reading {sample} from {file}.")
        df = pd.read_csv(file, sep='\t', compression='gzip').loc[:, [
            'Pos', 'ExonPos', 'Coverage']].rename(columns={'Coverage': sample})
        cov_df = cov_df.merge(df, on=['ExonPos', 'Pos'], how='outer')
    cov_df = cov_df.fillna(0).sort_values('ExonPos')
    cov_df['Chr'] = chrom
    # reorder columns
    cols = ['Chr', 'Pos', 'ExonPos'] + list(cov_df.columns)[2:-1]
    return cov_df.loc[:, cols]


def normalize_coverage(cov_df, norm_cov=100):
    norm_df = cov_df.set_index(['Chr', 'Pos', 'ExonPos'])
    norm_df = norm_df / norm_df.mean() * norm_cov
    return norm_df.reset_index()


def add_mean(norm_df):
    norm_df = norm_df.set_index(['Chr', 'Pos', 'ExonPos'])
    norm_df['meanCov'] = norm_df.mean(axis=1)
    norm_df['medianCov'] = norm_df.median(axis=1)
    norm_df['std'] = norm_df.std(axis=1)
    return norm_df.reset_index()


def filter_coverage(df, mincov=20, max_mean_std=20):
    filter_df = df.query('meanCov > @mincov and std < @max_mean_std')
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
    'sample_PON_path': '.',  # the folder with the sample/chrom coverages
    'output_path': '.',       # where the coverage files are saved
    'normCov': 100,       # to what value are coverages normalized#
    'minCov': 20,        # only exonPositions with the average coverage above minCov are kept
    'max_mean_std': 20,  # only exonPositions with a coverage std below max_mean_std are kept
    # only exonPositions straighing within std_factor * std around meanCoverage are kept
    'std_factor': 3,
}):
    '''
    load the coverages for all the PON files for one chromosome and write that to a file if given
    '''
    # load all sample coverages for one chromosome
    cov_df = combine_coverage(config['sample_PON_path'], chrom, sample_list)

    # normalize and add mean values and std
    mean_df = add_mean(normalize_coverage(cov_df, norm_cov=config['normCov']))

    # filter hard regions and outlying data points
    filter_df = filter_coverage(
        mean_df, mincov=config['minCov'], max_mean_std=config['max_mean_std'])
    final_df = remove_outliers(filter_df, std_factor=config['std_factor'])

    # output only if outpath is given
    if config['output_path']:
        file_name = f"{config['output_path']}/PON_coverage.{chrom}.removed.csv"
        print(f"Saving filtered PON file {file_name}")
        final_df.to_csv(file_name, sep='\t', index=False)
    return final_df
