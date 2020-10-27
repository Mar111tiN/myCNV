import pandas as pd
import numpy as np
from script_utils import show_output


def interpolate(df, data_col, ref_col='FullExonPos', expand_limit=20):
    cols = list(df.columns)
    # set FullExonPos as index for the interpolation method to work on proper intervals
    df = df.set_index(ref_col)
    df.loc[:, data_col] = df[data_col].interpolate(
        method='values', limit=expand_limit, limit_direction='both')
    return df.reset_index()[cols]


def one_col_rolling(df, df_filter, col, agg, window_size=200, expand_limit=20, normalize=False, debug=False):

    #
    org_cols = list(df.columns)
    # rolling left
    # get the right computation by passing agg to .agg()
    # only this allows passing methods as string
    df.loc[:, 'L'] = df_filter[col].rolling(window_size).agg(agg)
    # rolling right by shifting the L column
    df.loc[:, 'R'] = df.shift(-window_size + 1)['L']

    col_name = col + agg
    diff_name = col_name + "Diff"
    new_cols = org_cols + [col_name, diff_name]
    if debug:
        new_cols += [f'{col_name}L', f'{col_name}R']
    # skips interpolation if value == 0
    if interpolate:
        # interpolate missing values
        for c in ['L', 'R']:
            df = interpolate(df, c, expand_limit=expand_limit)

    # normalize values
    # not good for coverage and VAF
    if normalize:
        # normalize the data
        print('Normalizing data')
        _min = df['L'].min()
        _max = df['L'].max()
        for c in ['L', 'R']:
            df.loc[:, c] = (df[c] - _min) / (_max - _min)

    # get the Diff
    df.loc[:, diff_name] = ((df['R'] - df['L']) ** 2)
    df.loc[:, diff_name] = df[diff_name] / df[diff_name].max()
    # here, contribution of L and R is controlled by diff value
    df.loc[:, col_name] = df['R'] * \
        df[diff_name] + df['L'] * (1 - df[diff_name])

    # reduce to the right columns
    df = df.rename(columns=dict(L=f'{col_name}L', R=f'{col_name}R'))
    return df[new_cols]


def rolling_coverage(cov_df, config):
    '''
    cycle through the chroms and perform rolling window computations of data set in config
    '''

    # split the params dict for easier access
    params = config['coverage']
    filter_params = params['filter']
    data_params = params['data']
    chrom_dfs = []
    for chrom in cov_df['Chr'].unique():
        # restrict to chrom
        chrom_df = cov_df.query('Chr == @chrom').sort_values('FullExonPos')
        # get the params for filtering
        min_cov = filter_params['min_cov']
        min_PON_cov = filter_params['min_PON_cov']
        max_PON_std = filter_params['max_PON_std']

        # filter df
        filter_df = chrom_df.query(
            'Coverage >= @min_cov and PONmeanCov >= @min_PON_cov and PONstd < @max_PON_std')
        for data_col in data_params.keys():
            for agg in data_params[data_col].keys():
                window_size = data_params[data_col][agg]
                expand_limit = int(params['expand'] * window_size)
                # print(f"Computing rolling window for {agg} of {data_col} with window size {window_size} on {chrom}")
                chrom_df = one_col_rolling(chrom_df, filter_df, data_col, agg, window_size=window_size,
                                           expand_limit=expand_limit, normalize=params['normalize'], debug=config['debug'])
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')

    return df


def interpolate_fullexonpon(merge_df):
    chrom_dfs = []
    for chrom in merge_df['Chr'].unique():
        chrom_df = merge_df.query('Chr == @chrom')
        chrom_df = interpolate(chrom_df, 'FullExonPos',
                               ref_col='Pos', expand_limit=1000000)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')
    df.loc[:, 'FullExonPos'] = df['FullExonPos'].astype(int)
    return df


def mergeSNPnCov(cov_df, snp_df):

    # reduce the data to important columns
    # snp
    snp_keep_cols = list(snp_df.columns)[:3] + ['Depth', 'VAF']
    snp_df = snp_df.loc[:, snp_keep_cols]
    # cov
    cov_keep_cols = list(cov_df.columns)[
        :4] + ['log2ratiomean', 'log2ratiomeanDiff']
    cov_df = cov_df.loc[:, cov_keep_cols]

    # merge the data
    merge_df = cov_df.merge(snp_df, on=list(snp_df.columns[:3]), how='outer')

    # interpolate FullExonPos
    merge_df = interpolate_fullexonpon(merge_df)

    # interpolate the data
    for col in [col for col in merge_df.columns if 'log2ratio' in col]:
        merge_df = interpolate(merge_df, col, expand_limit=100)
    # reduce to VAF values
    snpcov_df = merge_df.query('VAF == VAF')
    return snpcov_df, cov_df


def apply_rolling_coverage(snp_df, cov_df, config):
    '''
    master function for rolling coverage
    '''

    cov_df = rolling_coverage(cov_df, config)

    snpcov_df, rolling_cov_df = mergeSNPnCov(cov_df, snp_df)

    return snpcov_df, rolling_cov_df


##################################################
#################################################
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
