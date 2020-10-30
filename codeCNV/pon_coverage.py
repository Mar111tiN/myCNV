import os
import pandas as pd
import numpy as np
from script_utils import show_output


def gather_PONcoverage_chrom(chrom, sample_list, config={}):
    cov_df = pd.DataFrame(columns=['ExonPos', 'Pos'])
    for sample in sample_list:
        sample_file = os.path.join(
            config['sample_PON_path'], f"{sample}.{chrom}.bedCov")
        if not os.path.isfile(sample_file):
            show_output(f"no file: {sample_file}", color="warning")
            continue
        if config['verbose_output']:
            show_output(f"Reading {sample} from {sample_file}.")
        df = pd.read_csv(sample_file, sep='\t', compression='gzip').loc[:, [
            'Pos', 'ExonPos', 'Coverage']].rename(columns={'Coverage': sample})
        cov_df = cov_df.merge(df, on=['ExonPos', 'Pos'], how='outer')
    cov_df = cov_df.fillna(0).sort_values('ExonPos')
    cov_df['Chr'] = chrom
    # reorder columns
    cols = ['Chr', 'Pos', 'ExonPos'] + list(cov_df.columns)[2:-1]
    return cov_df.loc[:, cols]


def gather_PONcoverage(chrom_list=[], sample_list=[], config={}):
    '''
    combine the PONcoverage for all chromosomes
    '''

    cov_dfs = []
    for chrom in chrom_list:
        show_output(f"Collecting PON coverages for {chrom}")
        cov_df = gather_PONcoverage_chrom(chrom, sample_list, config)
        cov_dfs.append(cov_df)
    cov_df_full = pd.concat(cov_dfs).reset_index(drop=True)
    return cov_df_full


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


def get_full_exon_pos(df):
    '''
    adds the accumulated exonic position (over all chroms)
    '''

    # save the output columns
    cols = list(df.columns)
    df = df.reset_index(drop=True)
    # adds the last ExonPos of chrom to start of next chromosome
    df.loc[:, 'chromStep'] = df.shift(1)['ExonPos'].fillna(0).astype(int)
    df.loc[df['Chr'] == df.shift(1)['Chr'], 'chromStep'] = 0
    df['chromAccum'] = df['chromStep'].cumsum()
    df['FullExonPos'] = df['ExonPos'] + df['chromAccum']
    cols = cols[:2] + ['FullExonPos'] + cols[2:]
    return df[cols]


def remove_outliers(df, std_factor=2.5):
    '''
    cycle through all sample cols, remove outliers
    '''

    for col in list(df.columns)[3:-3]:
        df.loc[np.abs(df['meanCov'] - df[col]) / df['std']
               > std_factor, col] = np.nan
    return add_mean(df.iloc[:, :-3])


def make_PON_coverage(sample_list, chrom_list=[f"chr{chrom + 1}" for chrom in range(22)] + ['chrX'], config={
    # provide a different chrom_list if you don't want standard ['chr1', 'chr2'...]
    'normCov': 100,       # to what value are coverages normalized#
    # only exonPositions straighing within std_factor * std around meanCoverage are kept
    'stdFactor': 3,
    'sample_PON_path': '.',
    'verbose_output': False  # the path to the bed_cover_file
}):
    # load all sample coverages for one chromosome
    show_output(
        f"Loading all PON coverages from {config['sample_PON_path']} for exome-wide normalization", time=True)
    cov_df = gather_PONcoverage(
        chrom_list=chrom_list, sample_list=sample_list, config=config)

    # normalize and add mean values and std
    show_output(f"Normalizing PON coverages", time=True)
    mean_df = add_mean(normalize_coverage(cov_df, norm_cov=config['normCov']))
    # add full exon coords to normalized PON coverage
    full_df = get_full_exon_pos(mean_df)
    show_output(f"Removing outliers", time=True)
    filter_df = remove_outliers(mean_df, std_factor=config['stdFactor'])
    # remove sample columns and addd full exon coords to filtered PON coverage
    filter_df = get_full_exon_pos(
        filter_df.loc[:, ['Chr', 'Pos', 'ExonPos', 'meanCov', 'medianCov', 'std']])
    show_output(f"Finished filtering of PON coverages",
                time=True, color="success")
    return full_df, filter_df
