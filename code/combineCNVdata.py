# try with pd merge
import os
import numpy as np
import pandas as pd
from script_utils import show_output


chrom_list = [f"chr{chrom + 1}" for chrom in range(22)] + ['chrX']


def combine_SNPdata(sample, sample_cnv_path="", verbose=False):
    snp_dfs = []
    file_base = os.path.join(sample_cnv_path, sample)
    for chrom in chrom_list:
        # reading SNP
        file = f"{file_base}.{chrom}.snp"
        if not os.path.isfile(file):
            show_output(f'No file {file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading SNP VAF from {chrom} of sample {sample} from {file}.")
        snp_df = pd.read_csv(file, sep='\t')
        snp_df[['Alt', 'AltDepth']] = snp_df['Alt'].str.extract(
            r"([AGCT])([0-9]+)")

        # reading snpEB
        file = f"{file_base}.{chrom}.snpEB"
        if not os.path.isfile(file):
            show_output(f'No file {file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading SNP EBscores from {chrom} of sample {sample} from {file}.")
        snpEB_df = pd.read_csv(file, sep='\t').loc[:, [
            'Chr', 'Start', 'Ref', 'Alt', 'EBscore', 'PoN-Alt']]
        snp_df = snp_df.merge(snpEB_df, on=['Chr', 'Start', 'Ref', 'Alt'])

        snp_dfs.append(snp_df)
    snp_df = pd.concat(snp_dfs)
    return snp_df.loc[:, ["Chr", "Start", "ExonPos", "Ref", "Depth", "Alt", "VAF", "EBscore", "PoN-Alt"]]


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
    cols = cols[:2] + ['FullExonPos'] + cols[2:] + ['chromAccum']
    return df[cols]


def combine_Covdata(sample, sample_cnv_path="", PON_cnv_path="", verbose=False):

    cover_dfs = []

    for chrom in chrom_list:
        # reading sampleCoverage
        sample_cov_file = os.path.join(
            sample_cnv_path, f"{sample}.{chrom}.cov")
        if not os.path.isfile(sample_cov_file):
            show_output(f'No file {file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading coverage from {chrom} of sample {sample} from {sample_cov_file}.")
        cov_df = pd.read_csv(sample_cov_file, sep='\t', compression="gzip")

        # reading PONcoverage
        pon_cov_file = os.path.join(PON_cnv_path, f"{chrom}.filtered.csv.gz")
        if not os.path.isfile(pon_cov_file):
            show_output(f'No file {file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading PON coverage of {chrom} from file {pon_cov_file}.")
        pon_df = pd.read_csv(pon_cov_file, sep='\t', compression="gzip").loc[:, [
            'Chr', 'Pos', 'ExonPos', 'meanCov', 'medianCov', 'std']]
        # column rename
        trans_dict = {col: f"PON{col}" for col in pon_df.columns[3:]}
        pon_df = pon_df.rename(columns=trans_dict)
        # merge sample with PON coverage
        sample_df = cov_df.merge(pon_df, on=['Chr', 'Pos', 'ExonPos'], how="inner").loc[:, [
            'Chr', 'Pos', 'ExonPos', 'Coverage', 'PONmeanCov', 'PONmedianCov', 'PONstd']]

        cover_dfs.append(sample_df)
    cover_df = pd.concat(cover_dfs)

    # normalize the coverage after concating
    cover_df['Coverage'] = cover_df['Coverage'] / \
        cover_df['Coverage'].mean() * 100
    cover_df['log2ratio'] = np.log2(
        cover_df['Coverage'] / cover_df['PONmeanCov'])
    return get_full_exon_pos(cover_df)


def get_full_exon_pos_from_cov(snp_df, cov_df):
    snp_cols = list(snp_df.columns)
    snp_df = snp_df.merge(cov_df.loc[:, ['Chr', 'chromAccum']].groupby(
        'Chr').first().reset_index(), on='Chr')
    snp_df['FullExonPos'] = snp_df['ExonPos'] + snp_df['chromAccum']
    cols = snp_cols[:2] + ['FullExonPos'] + snp_cols[2:]
    return snp_df[cols], cov_df.drop(columns=['chromAccum'])


def approx_log2ratio(snp_df, cov_df):
    '''
    takes the coverage data and approximates the log2ratio for that SNP from adjacent cov data
    '''

    # merge snp_df and cov_df and rename required columns
    merge_df = snp_df.merge(cov_df, on=[
        'Chr',
        'FullExonPos'
    ], how='outer').sort_values('FullExonPos').reset_index(drop=True).drop(
        columns=['Pos'] + list(cov_df.columns[4:-1])
    )

    # store the fitting SNPs
    merged_df = merge_df.query('EBscore == EBscore and log2ratio == log2ratio').drop(
        columns='ExonPos_y').rename(columns={'ExonPos_x': 'ExonPos'})

    # go on with the SNPs with missing log2ratio
    merge_df = merge_df.query('EBscore != EBscore or log2ratio != log2ratio').rename(columns={
        'ExonPos_x': 'ExonPos',
        'ExonPos_y': 'ExonPosL',
        'log2ratio': 'log2ratioL'
    })

    merge_dfs = []
    snp_cols = list(snp_df.columns) + ['log2ratio']

    for chrom in merge_df['Chr'].unique():
        merge = merge_df.query('Chr == @chrom')
        merge['ExonPosR'] = merge['ExonPosL'].fillna(method="bfill")
        merge['log2ratioR'] = merge['log2ratioL'].fillna(method="bfill")
        merge['ExonPosL'] = merge['ExonPosL'].fillna(method="ffill")
        merge['log2ratioL'] = merge['log2ratioL'].fillna(method="ffill")
        merge['log2ratio'] = merge['log2ratioL'] + (merge['log2ratioR'] - merge['log2ratioL']) / (
            merge['ExonPosR'] - merge['ExonPosL']) * (merge['ExonPos'] - merge['ExonPosL'])
        merge_dfs.append(merge.loc[:, snp_cols])
    snp_df = pd.concat(merge_dfs).sort_values(
        'FullExonPos').query('EBscore == EBscore')
    snp_df['log2ratio'] = snp_df['log2ratio'].fillna(method='bfill')
    return snp_df


def centerVAF(snp_df):
    '''
    attempting to correct for off-center VAF means
    '''

    # get the VAF mean
    meanVAF = snp_df.query('0.05 < VAF < 0.95')['VAF'].mean()
    # store the original VAF in orgVAF
    snp_df['orgVAF'] = snp_df['VAF']
    snp_df.loc[snp_df['VAF'] <= meanVAF,
               'VAF'] = snp_df['VAF'] / meanVAF * 0.5
    snp_df.loc[snp_df['VAF'] > meanVAF, 'VAF'] = 0.5 + \
        0.5 * (snp_df['VAF'] - meanVAF) / (1-meanVAF)
    return snp_df, meanVAF

# combine SNP data and covData


def get_covNsnp(sample, sample_cnv_path='', PON_cnv_path='', verbose=False, centerSNP=False):
    '''
    load the coverage_data for a sample and the heteroSNP data and apply the same fullExonCoords
    '''
    show_output(f'Loading coverage data for sample {sample}')
    cov_df = combine_Covdata(
        sample, sample_cnv_path=sample_cnv_path, PON_cnv_path=PON_cnv_path, verbose=verbose)
    show_output(f'Loading SNP data for sample {sample}')
    snp_df = combine_SNPdata(
        sample, sample_cnv_path=sample_cnv_path, verbose=verbose)
    # get full exonPos from cov_df and remove the fullAc
    snp_df, cov_df = get_full_exon_pos_from_cov(snp_df, cov_df)
    # # get lo
    # snp_df = approx_log2ratio(snp_df, cov_df)
    if centerSNP:
        snp_df, meanVAF = centerVAF(snp_df)
        show_output(
            f'Found SNPs offCenter at {meanVAF}. Centering SNPs around 0.5 {sample}')
    show_output(f"Finished loading sample {sample}", color="success")
    return snp_df, cov_df
