# try with pd merge
import os
import numpy as np
import pandas as pd
from script_utils import show_output
from codeCNV.rollingCov import apply_rolling_coverage
from codeCNV.rollingSNP import apply_rolling_SNP

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
        file = f"{file_base}.{chrom}.snp"
        if not os.path.isfile(file):
            show_output(f'No file {file}', color="warning")
            continue
        snp_dfs.append(snp_df)
    snp_df = pd.concat(snp_dfs).rename({'Start': 'Pos'}, axis=1)
    return snp_df.loc[:, ["Chr", "Pos", "ExonPos", "Ref", "Depth", "Alt", "VAF"]]


def combine_Covdata(sample, sample_cnv_path="", PON_cnv_path="", verbose=False, filtered=True):

    cover_dfs = []
    for chrom in chrom_list:
        # reading sampleCoverage
        sample_cov_file = os.path.join(
            sample_cnv_path, f"{sample}.{chrom}.cov")
        if not os.path.isfile(sample_cov_file):
            show_output(f'No file {sample_cov_file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading coverage from {chrom} of sample {sample} from {sample_cov_file}.")
        cov_df = pd.read_csv(sample_cov_file, sep='\t', compression="gzip")

        # reading PONcoverage
        full_or_filtered = "filtered" if filtered else "full"
        pon_cov_file = os.path.join(
            PON_cnv_path, f"{chrom}.{full_or_filtered}.csv.gz")
        if not os.path.isfile(pon_cov_file):
            show_output(f'No file {pon_cov_file}', color="warning")
            continue
        if verbose:
            show_output(
                f"Reading PON coverage of {chrom} from file {pon_cov_file}.")
        pon_df = pd.read_csv(pon_cov_file, sep='\t', compression="gzip").loc[:, [
            'Chr', 'Pos', 'FullExonPos', 'ExonPos', 'meanCov', 'medianCov', 'std']]
        # column rename PON columns to PON<col>
        trans_dict = {col: f"PON{col}" for col in pon_df.columns[4:]}
        pon_df = pon_df.rename(columns=trans_dict)
        # merge sample with PON coverage
        sample_df = cov_df.merge(pon_df, on=['Chr', 'Pos', 'ExonPos'], how="right").loc[:, [
            'Chr', 'Pos', 'FullExonPos', 'ExonPos', 'Coverage', 'PONmeanCov', 'PONmedianCov', 'PONstd']]

        # here recover missing FullExonPos from margin
        # get the off
        exon_start, full_start = sample_df.iloc[0][['ExonPos', 'FullExonPos']]
        offset = full_start - exon_start
        sample_df.loc[sample_df['FullExonPos'] != sample_df['FullExonPos'],
                      'FullExonPos'] = sample_df['ExonPos'] + offset
        sample_df.loc[:, 'FullExonPos'] = sample_df.loc[:,
                                                        'FullExonPos'].astype(int)
        cover_dfs.append(sample_df)
    # combine chrom data
    cover_df = pd.concat(cover_dfs)

    # normalize the coverage over the entire exome!
    cover_df['Coverage'] = cover_df['Coverage'].fillna(0)

    # take the mean without chrX because of male chrom!
    mean_cov = cover_df.loc[cover_df['Chr'] != "chrX", 'Coverage'].mean()
    cover_df.loc[:, 'Coverage'] = (cover_df['Coverage'] / mean_cov * 100)

    # loggable are the coverages, where log2ratio can be computed
    loggable = (cover_df['PONmeanCov'] * cover_df['Coverage'] != 0)
    cover_df.loc[loggable, 'log2ratio'] = np.log2(
        cover_df.loc[loggable, 'Coverage'] / cover_df.loc[loggable, 'PONmeanCov'])
    # mark regions without PON coverage as NAN
    cover_df.loc[~loggable, 'log2ratio'] = np.nan
    return cover_df

# combine SNP data and covData


def get_covNsnp(sample, sample_cnv_path='', PON_cnv_path='', verbose=False):
    '''
    load the coverage_data for a sample and the heteroSNP data and apply the same fullExonCoords
    '''
    show_output(f'Loading coverage data for sample {sample}')
    cov_df = combine_Covdata(
        sample, sample_cnv_path=sample_cnv_path, PON_cnv_path=PON_cnv_path, verbose=verbose)
    show_output(f'Loading SNP data for sample {sample}')
    snp_df = combine_SNPdata(
        sample, sample_cnv_path=sample_cnv_path, verbose=verbose)
    show_output(f"Finished loading sample {sample}", color="success")
    return cov_df, snp_df


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
    rolling_snpcov_df, cluster_df = apply_rolling_SNP(snpcov_df, config)
    show_output(f"Finished computations for sample {sample}.")
    if config['na_remove']:
        rolling_cov_df = rolling_cov_df.dropna()
        rolling_snpcov_df = rolling_snpcov_df.dropna()
        show_output(
            f"Removed missing values from rolling data for sample {sample}.")
    return cov_df, snp_df, rolling_cov_df, rolling_snpcov_df, cluster_df
