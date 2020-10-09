import os
import numpy as np
import pandas as pd
from script_utils import show_output


def get_approx_col_list(config):
    '''
    generates a list of dictionaries with data columns to be approximated
    from the config to be consumed by the approximator
    '''

    # approx_cols is the list of
    approx_cols = []
    snp_conf = config['heteroSNP']['windows']

    for col in snp_conf.keys():
        for mode in snp_conf[col].keys():
            approx_cols.append({
                'col': f"{col}{mode}",
                'pos_col': 'PosSNP',
                'trans_pos_col': 'PosCov'
            })
    cov_conf = config['coverage']['windows']
    for col in cov_conf.keys():
        for mode in cov_conf[col].keys():
            approx_cols.append({
                'col': f"{col}{mode}",
                'pos_col': 'PosCov',
                'trans_pos_col': 'PosSNP'
            })
    return approx_cols


def approximate_data(merge, col='VAFstd', pos_col='PosSNP', trans_pos_col='PosCov'):
    '''
    takes the data values from col at positions in pos_col 
    and linearly approximates data values into merged rows at positions in trans_pos_col
    '''

    cols = list(merge.columns)
    # find the adjacent positions for missing rows and store in PosL and PosR
    merge.loc[:, 'PosL'] = merge[pos_col].fillna(method="ffill")
    merge.loc[:, 'PosR'] = merge[pos_col].fillna(method="bfill")
    # find the adjacent data values for missing rows and store in L and R
    merge.loc[:, 'L'] = merge[col].fillna(method="ffill")
    merge.loc[:, 'R'] = merge[col].fillna(method="bfill")
    # approximate the missing values
    merge.loc[merge[col] != merge[col], col] = merge['L'] + (merge['R'] - merge['L']) / (
        merge['PosR'] - merge['PosL']) * (merge[trans_pos_col] - merge['PosL'])
    # close the gaps
    merge.loc[:, col] = merge[col].fillna(
        method='bfill').fillna(method='ffill')

    # return the only the original columns with filled in values
    return merge.loc[:, cols]


def mergeSNP2Cov(snp_df, cov_df, config):
    '''
    for clustering, all data points from SNP
    '''
    # get the required snp_df columns for merge
    snp_cols = ['Chr', 'Start', 'FullExonPos', 'ExonPos', 'VAF']
    for metrix in ['absVAFsum', 'deltaVAFstd', 'VAFstd']:
        snp_cols += [col for col in snp_df.columns if col.startswith(metrix)]

    snp_chrom_df = snp_df.loc[:, snp_cols]

    # get the required cov_df columns for merge
    cov_chrom_df = cov_df.loc[:, [
        'Chr', 'Pos', 'FullExonPos', 'ExonPos', 'log2ratiomeanDiff', 'log2ratiomean']]

    # do the merge and rename the respective ExonPos
    merge_df = snp_chrom_df.merge(cov_chrom_df, on=['Chr', 'FullExonPos'], how='outer').sort_values('FullExonPos').rename(columns={
        'ExonPos_x': 'PosSNP',
        'ExonPos_y': 'PosCov'
    })

    # merge chromosomal start coords
    merge_df.loc[merge_df['Start'] !=
                 merge_df['Start'], 'Start'] = merge_df['Pos']
    merge_df.loc[:, 'Start'] = merge_df['Start'].astype(int)
    merge_df = merge_df.drop(columns='Pos').reset_index(
        drop=True).sort_values('FullExonPos')

    # store the fitting values as merged_df
    merged_df = merge_df.query('VAF == VAF and log2ratiomean == log2ratiomean')

    # go on with the SNPs with non-fitting data
    merge_df = merge_df.query('VAF != VAF or log2ratiomean != log2ratiomean')

    # get the data columns for the approximator
    data_col_list = get_approx_col_list(config)

    # go through the chromosomes and do the approximation
    merge_dfs = []
    for chrom in merge_df['Chr'].unique():
        chrom_merge_df = merge_df.query('Chr == @chrom')
        for data in data_col_list:
            chrom_merge_df = approximate_data(
                chrom_merge_df, col=data['col'], pos_col=data['pos_col'], trans_pos_col=data['trans_pos_col'])
            chrom_merge_df = approximate_data(
                chrom_merge_df, col=data['col']+"Diff", pos_col=data['pos_col'], trans_pos_col=data['trans_pos_col'])
        merge_dfs.append(chrom_merge_df)
    # concat the chroms and add the already merged df
    merge_df = pd.concat(merge_dfs + [merged_df]).sort_values(
        'FullExonPos').rename(columns={'PosSNP': 'ExonPos'})
    # transfer the missing positions from
    merge_df.loc[:, 'ExonPos'] = merge_df['ExonPos'].fillna(merge_df['PosCov'])
    return merge_df.drop(columns='PosCov')
