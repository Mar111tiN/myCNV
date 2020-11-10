import pandas as pd
import numpy as np
from script_utils import show_output
from codeCNV.combineCNVdata import get_covNsnp


def interpolate(df, data_col, ref_col='FullExonPos', expand_limit=20):
    '''
    interpolates missing values in data_col using linear interpolation based on ref_col
    '''

    cols = list(df.columns)
    # set FullExonPos as index for the interpolation method to work on proper intervals
    df = df.reset_index(drop=False).set_index(ref_col, drop=False)
    df.loc[:, data_col] = df[data_col].interpolate(
        method='values', limit=expand_limit, limit_direction='both')
    return df.set_index('index')[cols]


def one_col_rolling(df, df_filter, col, aggr, window_size=200, expand_limit=20, normalize=False, debug=False, diff_exp=2, ddof=0):
    '''
    performs rolling computation of <agg> on data column <col> with given window size
    the aggregation has to be a string expression understood by the agg-function of the pandas.groupby API
    computation is performed on a left and right rolling window
    missing margins are filled by the counterpart window function
    a diff column is included ()

    '''

    org_cols = list(df.columns)
    # rolling left
    # get the right computation by passing aggr to .agg()
    # only this allows passing methods as string
    df.loc[:, 'L'] = df_filter[col].rolling(window_size).agg(aggr, ddof=ddof)
    # rolling right by shifting the L column
    df.loc[:, 'R'] = df.shift(-window_size + 1)['L']

    col_name = col + aggr
    diff_name = col_name + "Diff"
    new_cols = org_cols + [col_name, diff_name]
    if debug:
        new_cols += [f'{col_name}L', f'{col_name}R']
    # skips interpolation if value == 0
    if interpolate:
        # interpolate missing values
        for c in ['L', 'R']:
            df = interpolate(df, c, expand_limit=expand_limit)
    # fill the margins
    L_margin = df['L'].first_valid_index()
    df.loc[:L_margin, 'L'] = df['R']
    R_margin = df['R'].last_valid_index() + 1
    df.loc[R_margin:, 'R'] = df['L']

    # get the Diff
    df.loc[:, diff_name] = np.abs(df['R'] - df['L'])
    # normalize to max
    df.loc[:, diff_name] = df[diff_name] / df[diff_name].max()
    # here, contribution of L and R is controlled by diff value
    df.loc[:, col_name] = df['R'] * \
        df[diff_name] + df['L'] * (1 - df[diff_name])

    # square the diff
    df.loc[:, diff_name] = df[diff_name] ** diff_exp

    if debug:
        # specify col names of L and R
        df = df.rename(columns=dict(L=f'{col_name}L', R=f'{col_name}R'))

    # reduce to the right columns
    return df[new_cols]


def llh(data, mean, sigma):
    '''
    compute the density function for a given gaussian
    takes a pd.Series or np.array
    '''
    s = np.sqrt(2 * np.pi) * sigma
    return np.exp((data - mean)**2 / (-2*(sigma**2))) / s


def compute_coverage_llh(df, config):
    '''
    computes the local log-likelihood of belonging to the center gaussian
    '''

    # get config params
    params = config['coverage']['llh']

    min_log2ratio, max_log2ratio = params['center_range']
    # get the sigma and mean of the center band log2ratio
    center_logs = df.query('@min_log2ratio < log2ratio < @max_log2ratio')[
        'log2ratio']
    sigma = center_logs.std() * params['sigma_factor']
    mean = center_logs.mean()
    show_output(
        f"Computing log-likelihood of log2ratio belonging to center gaussian [mean:{round(mean, 3)}, sigma:{round(sigma,3)}]")
    df.loc[:, 'covllh'] = llh(df['log2ratio'], mean, sigma)

    return df


def rolling_coverage(cov_df, config):
    '''
    cycle through the chroms and perform rolling window computations of data set in config
    '''

    # split the params dict for easier access
    params = config['coverage']
    filter_params = params['filter']
    data_params = params['rolling_data']
    debug = config['debug']
    # get the params for filtering
    min_cov = filter_params['min_cov']
    min_PON_cov = filter_params['min_PON_cov']
    max_PON_std = filter_params['max_PON_std']
    chrom_dfs = []
    for chrom in cov_df['Chr'].unique():
        # restrict to chrom
        chrom_df = cov_df.query('Chr == @chrom').sort_values('FullExonPos')
        # filter df
        filter_df = chrom_df.query(
            'Coverage >= @min_cov and PONmeanCov >= @min_PON_cov and PONstd < @max_PON_std')
        for data_col in data_params.keys():
            for agg in data_params[data_col].keys():
                window_size = data_params[data_col][agg]
                expand_limit = int(params['expand'] * window_size)
                # show_output(f"Computing rolling window for {agg} of {data_col} with window size {window_size} on {chrom}")
                chrom_df = one_col_rolling(chrom_df, filter_df, data_col, agg,
                                           window_size=window_size,
                                           expand_limit=expand_limit,
                                           normalize=params['normalize'],
                                           diff_exp=config['diff_exp'],
                                           debug=debug)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')

    # now do global normalization for sum aggregations:
    # cycle through rolling_data
    for data_col in data_params.keys():
        for agg in data_params[data_col].keys():
            # only do normalization for sum aggregations
            if not agg == "sum":
                continue
            show_output(f"Normalizing {data_col} {agg}")
            # get the columns for normalization
            col_name = data_col + agg
            cols = [col_name]
            if debug:
                cols += [f'{col_name}L', f'{col_name}R']

            for c in cols:

                _min = df[c].min()
                _max = df[c].max()
                df.loc[:, c] = (df[c] - _min) / (_max - _min)
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
    snp_keep_cols = list(snp_df.columns)[:3] + ['Depth', 'EBscore', 'VAF']
    snp_df = snp_df.loc[:, snp_keep_cols]
    # cov
    cov_keep_cols = list(cov_df.columns)[:4]
    for data in ['log2ratio', 'covllh']:
        cov_keep_cols += [col for col in cov_df.columns if data in col]

    cov_df = cov_df.loc[:, cov_keep_cols]

    # merge the data
    merge_df = cov_df.merge(snp_df, on=list(snp_df.columns[:3]), how='outer')

    # interpolate FullExonPos
    merge_df = interpolate_fullexonpon(merge_df)

    # interpolate the data
    for col in [col for col in merge_df.columns if 'log2ratio' in col or 'covllh' in col]:
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
    cov_df = compute_coverage_llh(cov_df, config)

    cov_df = rolling_coverage(cov_df, config)

    snpcov_df, rolling_cov_df = mergeSNPnCov(cov_df, snp_df)

    return snpcov_df, rolling_cov_df


def expand_SNPdata(snp_df, config):
    '''
    retrieve a few data columns locally to use rolling windows on
    this needs to be done chromosome-wise in order to avoid gap effects
    VAF limits are also applied here
    '''

    # split the params dict for easier access
    params = config['heteroSNP']
    filter_params = params['filter']
    data_params = params['data']
    # reduce the snp_df using config limits
    VAFmin, VAFmax = filter_params['VAF']
    snp_df = snp_df.query('@VAFmin < VAF < @VAFmax')

    # get the new features from VAFs
    snp_df.loc[:, 'absVAF'] = np.abs(snp_df['VAF'] - 0.5) * 2
    # get the local VAF difference chrom based
    dfs = []
    for chrom in snp_df['Chr'].unique():
        chrom_df = snp_df.query('Chr == @chrom')
        chrom_df.loc[:, 'deltaVAF'] = np.abs(
            chrom_df['VAF'] - chrom_df.shift(1)['VAF']).fillna(0)
        dfs.append(chrom_df)
    snp_df = pd.concat(dfs).sort_values('FullExonPos')
    return snp_df.reset_index(drop=True)


def rolling_SNP(snp_df, config):
    '''
    cycle through the chroms and perform rolling window computations of snp data set in config
    '''

    # split the params dict for easier access
    params = config['heteroSNP']
    filter_params = params['filter']
    data_params = params['rolling_data']
    # reduce the snp_df using config limits
    VAFmin, VAFmax = filter_params['VAF']
    minDepth = filter_params['minDepth']

    # cycle through chroms for
    chrom_dfs = []
    for chrom in snp_df['Chr'].unique():
        # restrict to chrom
        chrom_df = snp_df.query('Chr == @chrom').sort_values('FullExonPos')
        # filter df
        # .query('@VAFmin < VAF < @VAFmax and
        filter_df = snp_df.query(
            'Depth >= @minDepth')
        for data_col in data_params.keys():
            for agg in data_params[data_col].keys():
                window_size = data_params[data_col][agg]
                expand_limit = int(params['expand'] * window_size)
                # show_output(f"Computing rolling window for {agg} of {data_col} with window size {window_size} on {chrom}")
                chrom_df = one_col_rolling(chrom_df, filter_df, data_col, agg, window_size=window_size,
                                           expand_limit=expand_limit, normalize=params['normalize'], debug=config['debug'], ddof=config['ddof'])
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).sort_values('FullExonPos')
    return df


def apply_rolling_SNP(snp_df, config):
    '''
    expands the SNP data to absVAF and deltaVAF and performs rolling window computations
    set in the config
    '''

    # get extra data
    snp_df = expand_SNPdata(snp_df, config)
    # do the rolling
    snp_df = rolling_SNP(snp_df, config)
    return snp_df


def rollingCNV(sample, sample_cnv_path, PON_cnv_path, config):
    '''
    combines all the hetSNP and coverage data per sample and 
    performs rolling computations for clustering
    returns the combined raw data and the (optionally na_removed) rolling data 
    '''

    # combine the chromosome data and associate coverage data with pon coverage
    snp_df, cov_df = get_covNsnp(
        sample,
        sample_cnv_path=sample_cnv_path,
        PON_cnv_path=PON_cnv_path,
        verbose=config['debug']
    )

    # apply rolling coverage
    show_output(
        f"Performing rolling coverage computation for sample {sample}.")
    snpcov_df, rolling_cov_df = apply_rolling_coverage(snp_df, cov_df, config)

    # apply rolling SNP
    show_output(
        f"Performing rolling computation for hetSNP data of sample {sample}.")
    rolling_snpcov_df = apply_rolling_SNP(snpcov_df, config)
    show_output(f"Finished computations for sample {sample}.")
    if config['na_remove']:
        rolling_cov_df = rolling_cov_df.dropna()
        rolling_snpcov_df = rolling_snpcov_df.dropna()
        show_output(
            f"Removed missing values from rolling data for sample {sample}.")
    return cov_df, snp_df, rolling_cov_df, rolling_snpcov_df
