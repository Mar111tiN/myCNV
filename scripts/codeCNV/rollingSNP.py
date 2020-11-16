import pandas as pd
import numpy as np

from codeCNV.rollingCNV import llh, rolling_data, one_col_rolling, get_CNV_blocks, interpolate
from script_utils import show_output


def make_get_density(window_size=20):
    '''
    higher order functional helper for returning a density computer for given window_size
    '''

    def SNPdensity(data):
        return (data.max() - data.min()) / window_size
    return SNPdensity


def compute_snp_llh(df, mean=0.5, sigma=0.2):
    '''
    computes the local log-likelihood of belonging to the center gaussian
    '''

    show_output(
        f"Computing log-likelihood of VAF belonging to center gaussian [mean:{round(mean, 3)}, sigma:{round(sigma,3)}]")
    df.loc[:, 'snpLLH'] = llh(df['VAF'], mean, sigma)

    # for homoSNPs reduce the VAFs to the ones above mean
    upper_vafs = df.query('@mean < VAF')['VAF']
    # then compute the hsnpLLH

    show_output(
        f"Computing log-likelihood of VAF belonging to purity100  [mean:1, sigma:{round(sigma,3)}]")
    # these are called hsnp
    # upper_vafs only contains half the snps, the remaining have to be interpolated
    df.loc[:, 'hsnpLLH'] = llh(upper_vafs, 1, sigma)
    df = interpolate(df, 'hsnpLLH', expand_limit=50)
    return df


def remove_fallSNP(snp_df, mean=0.5, std=0.2, params={}):
    '''
    removes the falling SNP probably caused by mismapping
    '''

    window = params['offVAFwindow']
    cutoff = params['maxFallSNP']

    # get the density computer for rolling
    get_SNPdensity = make_get_density(window)
    # cycle through chroms
    chrom_dfs = []
    for chrom in snp_df['Chr'].unique():
        df = snp_df.query('Chr == @chrom').reset_index(drop=True)
        if len(df.index) < 10:
            continue
        # get the snp
        df = one_col_rolling(df, df.query(
            'VAF < 0.95'), 'ExonPos', get_SNPdensity, window_size=window, diff_exp=4)
        df.loc[:, 'SNPdensity'] = df['SNPdensity'] / df['SNPdensity'].mean()

        # get the offVAFsum
        df = one_col_rolling(df, df.query(
            'VAF < 0.95'), 'offVAF', 'sum', window_size=window, normalize=True, diff_exp=4)

        # combine both metrices
        df.loc[:, 'fallSNP'] = df['SNPdensity'] * df['offVAFsum']
        # now remove the ones below average VAFstd
        df = df.query('VAF > @mean - @std / 2 or fallSNP > @cutoff')
        chrom_dfs.append(df)

    return pd.concat(chrom_dfs).sort_values('FullExonPos').reset_index(drop=True)


def expand_SNPdata(snp_df, config):
    '''
    retrieve a few data columns locally to use rolling windows on
    this needs to be done chromosome-wise in order to avoid gap effects
    VAF limits are also applied here
    '''

    # split the params dict for easier access
    params = config['snp']
    filter_params = params['filter']
    # data_params = params['data']

    # reduce the snp_df using lower config limit
    # upper limit has to be set later as we still need the homoSNP llh
    VAFmin, VAFmax = filter_params['VAF']
    snp_df = snp_df.query('@VAFmin < VAF')

    # get std and mean of VAF
    minVAF, maxVAF = params['LLH']['center_range']
    # get the sigma and mean of the center band VAF (extracted as pd.Series center_vafs)
    center_vafs = snp_df.query('@minVAF < VAF < @maxVAF')['VAF']
    # get width of gaussian from std * sigma_factor
    VAFstd = center_vafs.std()
    VAFmean = center_vafs.mean()

    # get additional features from VAFs
    snp_df.loc[:, 'offVAF'] = (snp_df['VAF'] - VAFmean) * 2
    # absolute values for cluster
    snp_df.loc[:, 'absVAF'] = np.abs(snp_df['offVAF'])

    ########## remove fallSNP ########
    fs_params = params['fallSNP']
    if fs_params['run']:
        show_output('Removing falling SNPs')
        snp_df = remove_fallSNP(snp_df, mean=VAFmean,
                                std=VAFstd, params=fs_params)

    ######## LLH  #####################
    # get the snpLLH and hsnpLLH
    # get config params
    sigma = VAFstd * params['LLH']['sigma_factor']
    # hsnpLLH is computed in order to rescue high absVAF that would have been filtered out
    # lower VAF is already removed because density of VAF ~0 is highly irregular and would confound density estimates
    snp_df = compute_snp_llh(snp_df, mean=VAFmean, sigma=sigma)
    return snp_df.query('VAF < @VAFmax').reset_index(drop=True)


def rolling_SNP(snp_df, config):
    '''
    cycle through the chroms and perform rolling window computations of snp data set in config
    '''

    # split the params dict for easier access
    params = config['snp']
    filter_params = params['filter']
    data_params = params['rolling_data']
    debug = config['debug']

    minDepth = filter_params['minDepth']
    filter_df = snp_df.query('Depth >= @minDepth')
    show_output("Performing rollingSNP computations.")
    rolling_df = rolling_data(
        snp_df, filter_df, expand=params['expand'], ddof=config['ddof'], debug=debug, data_params=data_params)
    return rolling_df


def apply_rolling_SNP(snp_df, config):

    # get extra data
    snp_df = expand_SNPdata(snp_df, config)
    # do the rolling
    snp_df = rolling_SNP(snp_df, config)
    # get the CNV and Center blocks
    snp_df = get_CNV_blocks(snp_df, 'snpLLH', config)

    # select columns for output
    base_cols = list(snp_df.columns[:4])

    snp_cols = [col for col in snp_df.columns[4:]
                if not 'log2' in col and not 'cov' in col and not 'off' in col]
    rolling_snp_df = snp_df[base_cols + snp_cols]
    cluster_cols = ['log2ratio', 'log2ratiomean',
                    'VAF', 'absVAF', 'absVAFmean']
    cluster_cols += [col for col in snp_df.columns if 'Center' in col or 'CNV' in col]
    cluster_df = snp_df[base_cols + cluster_cols]
    return rolling_snp_df, cluster_df
