import pandas as pd
import numpy as np

from codeCNV.rollingCNV import llh, rolling_data, get_CNV_blocks, interpolate
from script_utils import show_output


def compute_coverage_llh(df, config):
    '''
    computes the local log-likelihood of belonging to the center gaussian
    '''

    # get config params
    params = config['cov']['LLH']

    min_log2ratio, max_log2ratio = params['center_range']
    # get the sigma and mean of the center band log2ratio
    center_logs = df.query('@min_log2ratio < log2ratio < @max_log2ratio')[
        'log2ratio']
    sigma = center_logs.std() * params['sigma_factor']
    mean = center_logs.mean()
    show_output(
        f"Computing log-likelihood of log2ratio belonging to center gaussian [mean:{round(mean, 3)}, sigma:{round(sigma,3)}]")
    df.loc[:, 'covLLH'] = llh(df['log2ratio'], mean, sigma)

    return df


def rolling_coverage(cov_df, config):
    '''
    cycle through the chroms and perform rolling window computations of data set in config
    '''

    # split the params dict for easier access
    params = config['cov']
    filter_params = params['filter']
    data_params = params['rolling_data']

    # get the params for filtering
    min_cov = filter_params['min_cov']
    min_PON_cov = filter_params['min_PON_cov']
    max_PON_std = filter_params['max_PON_std']

    cov_df = cov_df.sort_values('FullExonPos')
    filter_df = cov_df.query(
        'Coverage >= @min_cov and PONmeanCov >= @min_PON_cov and PONstd < @max_PON_std')

    cov_df = rolling_data(
        cov_df, filter_df, expand=params['expand'], ddof=config['ddof'], debug=config['debug'], data_params=data_params)

    return cov_df


def interpolate_fullexonpon(merge_df):
    chrom_dfs = []
    for chrom in merge_df['Chr'].unique():
        chrom_df = merge_df.query('Chr == @chrom')
        chrom_df = interpolate(chrom_df, 'FullExonPos',
                               ref_col='Pos', expand_limit=1000000)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')
    df.loc[:, 'FullExonPos'] = df['FullExonPos'].fillna(0).astype(int)
    return df


def mergeSNPnCov(cov_df, snp_df):

    # reduce the data to important columns
    # snp
    snp_keep_cols = list(snp_df.columns)[:3] + ['Depth', 'VAF']
    snp_df = snp_df.loc[:, snp_keep_cols]
    # cov
    snpcov_keep_cols = list(cov_df.columns)[:4]
    # columns for snpcov_df only need to contain the basics
    for data in ['log2ratio', 'covC']:
        snpcov_keep_cols += [col for col in cov_df.columns if data in col]
    # columns for cov_df should contain most data
    cov_keep_cols = snpcov_keep_cols + \
        [col for col in cov_df.columns if 'covLLH' in col]

    cov_df = cov_df.loc[:, cov_keep_cols]

    merge_df = cov_df.loc[:, snpcov_keep_cols]
    # merge the data
    merge_df = merge_df.merge(snp_df, on=list(
        snp_df.columns[:3]), how='outer')

    # interpolate FullExonPos
    merge_df = interpolate_fullexonpon(merge_df)
    # interpolate the data for all added fields
    for col in [col for col in merge_df.columns if 'log2ratio' in col or 'covC' in col]:
        merge_df = interpolate(merge_df, col, expand_limit=100)
    # reduce to VAF values
    snpcov_df = merge_df.query('VAF == VAF')

    cov_df = cov_df.query('log2ratiomean == log2ratiomean')
    return snpcov_df, cov_df


def apply_rolling_coverage(snp_df, cov_df, config):
    '''
    master function for rolling coverage
    '''
    # reduce cov_df to valid data
    cov_df = cov_df.query('log2ratio == log2ratio')

    # compute llh
    show_output(
        f"Computing covCenter log-likelihood.")
    cov_df = compute_coverage_llh(cov_df, config)
    # compute llh
    show_output(
        f"Performing rolling coverage.")
    cov_df = rolling_coverage(cov_df, config)
    show_output(
        f"Identifying CNV blocks.")
    cov_df = get_CNV_blocks(cov_df, 'covLLH', config)
    snpcov_df, rolling_cov_df = mergeSNPnCov(cov_df, snp_df)
    return snpcov_df, rolling_cov_df
