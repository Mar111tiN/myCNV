import pandas as pd
import numpy as np
from script_utils import show_output


def expand_SNPdata(snp_df, config):
    '''
    retrieve a few data columns locally to use rolling windows on
    this needs to be done chromosome-wise in order to avoid gap effects
    '''

    # reduce the snp_df using config limits
    VAFmin = config['VAFlimits'][0]
    VAFmax = config['VAFlimits'][1]
    minDepth = config['minDepth']
    minEBscore = config['minEBscore']
    df = snp_df.query(
        '@VAFmin < VAF < @VAFmax and Depth >= @minDepth and EBscore > @minEBscore')

    # find the center
    center = df['VAF'].mean()
    show_output(f"heteroSNP centered around {center}", time=False)
    # offCenter
    df.loc[:, 'absVAF'] = np.abs(df['VAF'] - center)

    # get the local VAF difference chrom based
    dfs = []
    for chrom in df['Chr'].unique():
        chrom_df = df.query('Chr == @chrom')
        chrom_df['deltaVAF'] = np.abs(
            chrom_df['VAF'] - chrom_df.shift(1)['VAF']).fillna(0)
        dfs.append(chrom_df)
    snp_df = pd.concat(dfs).sort_values('FullExonPos')

    return snp_df


def get_cols(col, agg='mean', modes=['L', 'R', 'Diff', '']):
    '''
    creates for each col a dict for looped computation
    {'L': 'VAVsumL', 'R': 'VAVsumR', 'Diff': 'VAVsumDiff', '': 'VAFsum'}}
    '''
    cols = {mode: col + agg + mode for mode in modes}
    return cols


def get_rolling_metrix_chrom(df, col='VAF', agg='sum', chrom='', window_size=20):
    '''
    take a column and an aggregation produce rolling windows from it for each chromosome
    '''

    df = df.query('Chr == @chrom')
    cols = get_cols(col, agg)

    # get the right computation
    if agg == 'std':
        df.loc[:, cols['L']] = df.rolling(window_size)[col].std()
    if agg == 'sum':
        df.loc[:, cols['L']] = df.rolling(window_size)[col].sum()
    if agg == 'mean':
        df.loc[:, cols['L']] = df.rolling(window_size)[col].mean()

    # get the right window by shifting the left
    df.loc[:, cols['R']] = df.shift(-window_size + 1)[cols['L']]
    # fillup the margins
    df.loc[:, cols['L']] = df[cols['L']].fillna(method='bfill')
    df.loc[:, cols['R']] = df[cols['R']].fillna(method='ffill')
    return df


def get_rolling_metrix(df, col='VAF', agg='mean', window_size=20, normalize=True, keep_LR=False):
    '''
    wrapper to apply get_rolling_metrix_chrom per chromosome
    '''
    org_cols = list(df.columns)
    chrom_dfs = []
    for chrom in df['Chr'].unique():
        chrom_df = get_rolling_metrix_chrom(
            df, col=col, agg=agg, chrom=chrom, window_size=window_size)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')

    cols = get_cols(col, agg)
    if normalize:
        # normalize the data
        show_output('Normalizing data', time=False)
        _min = df[cols['L']].min()
        _max = df[cols['L']].max()
        for side in ['L', 'R']:
            c = cols[side]
            df[c] = (df[c] - _min) / (_max - _min)
    # get the Diff
    df[cols['Diff']] = ((df[cols['L']] - df[cols['R']]) / 2) + 0.5
    df[cols['']] = df[cols['L']] * df[cols['Diff']] + \
        df[cols['R']] * (1 - df[cols['Diff']])
    # reduce the columns

    keep_cols = org_cols + [cols['Diff'], cols['']]
    if keep_LR:
        keep_cols += [cols['L'], cols['R']]

    return df[keep_cols]


def rolling_it(df, config):

    windows = config['windows']
    for col in windows.keys():
        for agg in windows[col].keys():
            window_size = windows[col][agg]
            show_output(
                f"Computing rolling window for {agg} of {col} with window size {window_size}", time=False)
            df = get_rolling_metrix(
                df, col=col, agg=agg, window_size=window_size, normalize=config['normalize'])
    return df


def add_rolling_data(snp_df, cov_df, config):
    '''
    add the rolling metrices needed to get local data
    '''

    # add extra cols to snp_df
    snp_df = expand_SNPdata(snp_df, config['heteroSNP'])

    # get the rolling metrices for snp_df
    snp_df = rolling_it(snp_df, config['heteroSNP'])
    # get the rolling metrices for cov_df
    cov_df = rolling_it(cov_df, config['coverage'])

    return snp_df, cov_df
