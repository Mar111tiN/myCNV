# try with pd merge
import os
import re
import numpy as np
import pandas as pd
from script_utils_CNV import show_output


def combine_coverage(cov_list, pon_cov_path="", verbose=False):
    """
    reading coverage data for all chroms and normalize with respective PON coverage
    """

    # get the sample for verbose output
    sample = os.path.basename(cov_list[0].split(".")[0])
    # chrom_pat for extracting the chrom pattern from the file list
    chrom_pat = re.compile(r"\.(chr[0-9X]+)\.")

    cover_dfs = []
    for sample_cov_file in cov_list:
        # extracting chrom from sample_cov_file
        chrom = re.search(chrom_pat, sample_cov_file).group(1)

        # reading sampleCoverage
        if not os.path.isfile(sample_cov_file):
            show_output(f"No file {sample_cov_file}", color="warning")
            continue
        if verbose:
            show_output(
                f"Reading coverage from {chrom} of sample {sample} from {sample_cov_file}."
            )
        cov_df = pd.read_csv(sample_cov_file, sep="\t", compression="gzip")
        # reading PONcoverage
        #!!! ExonPos seems not to be compatible anymore
        pon_cov_file = os.path.join(pon_cov_path, f"{chrom}.filtered.csv.gz")
        if not os.path.isfile(pon_cov_file):
            show_output(f"No file {pon_cov_file}", color="warning")
            continue
        if verbose:
            show_output(f"Reading PON coverage of {chrom} from file {pon_cov_file}.")

        ###############################################
        # reading PONcoverage
        #!!! ExonPos seems not to be compatible anymore
        # --> ExonPos removed
        pon_df = pd.read_csv(pon_cov_file, sep="\t", compression="gzip").loc[
            :, ["Chr", "Pos", "FullExonPos", "meanCov", "medianCov", "std"]
        ]
        # column rename PON columns to PON<col>
        trans_dict = {col: f"PON{col}" for col in pon_df.columns[3:]}
        pon_df = pon_df.rename(columns=trans_dict)
        # merge sample with PON coverage
        ###############################################
        # reading PONcoverage
        #!!! ExonPos seems not to be compatible anymore
        sample_df = cov_df.merge(pon_df, on=["Chr", "Pos"], how="right").loc[
            :,
            [
                "Chr",
                "Pos",
                "FullExonPos",
                "ExonPos",
                "Coverage",
                "PONmeanCov",
                "PONmedianCov",
                "PONstd",
            ],
        ]
        # here recover missing FullExonPos from margin
        # get the offset from ExonPos relative to FullExonPos to infer missing positions
        # ########!!!!!######### check if this is working good!!!
        exon_start, full_start = sample_df.iloc[0][["ExonPos", "FullExonPos"]]
        offset = full_start - exon_start
        sample_df.loc[
            sample_df["FullExonPos"] != sample_df["FullExonPos"], "FullExonPos"
        ] = (sample_df["ExonPos"] + offset)
        sample_df.loc[:, "FullExonPos"] = sample_df.loc[:, "FullExonPos"].astype(int)
        cover_dfs.append(sample_df)
    # combine chrom data
    cover_df = pd.concat(cover_dfs)

    # normalize the coverage over the entire exome!
    # ########!!!!!######### maybe sometimes there are too many 0 which is weird for normalization
    cover_df["Coverage"] = cover_df["Coverage"].fillna(0)

    # take the mean without chrX because of male chrom!
    mean_cov = cover_df.loc[cover_df["Chr"] != "chrX", "Coverage"].mean()
    cover_df.loc[:, "Coverage"] = cover_df["Coverage"] / mean_cov * 100

    # loggable are the coverages, where log2ratio can be computed
    loggable = cover_df["PONmeanCov"] * cover_df["Coverage"] != 0
    cover_df.loc[loggable, "log2ratio"] = np.log2(
        cover_df.loc[loggable, "Coverage"] / cover_df.loc[loggable, "PONmeanCov"]
    )
    # mark regions without PON coverage as NAN
    cover_df.loc[~loggable, "log2ratio"] = np.nan
    return cover_df


def combine_SNP(snp_list, verbose=False):
    """
    reading hetSNP data for all chroms
    """

    # get the sample for verbose output
    sample = os.path.basename(snp_list[0].split(".")[0])
    # chrom_pat for extracting the chrom pattern from the file list
    chrom_pat = re.compile(r"\.(chr[0-9X]+)\.")
    snp_dfs = []
    for snp_file in snp_list:
        # extracting chrom from sample_cov_file
        chrom = re.search(chrom_pat, snp_file).group(1)
        # reading SNP
        if not os.path.isfile(snp_file):
            show_output(f"No file {snp_file}", color="warning")
            continue
        if verbose:
            show_output(
                f"Reading SNP VAF from {chrom} of sample {sample} from {snp_file}."
            )
        snp_df = pd.read_csv(snp_file, sep="\t")
        snp_df[["Alt", "AltDepth"]] = snp_df["Alt"].str.extract(r"([AGCT])([0-9]+)")

        snp_dfs.append(snp_df)
    snp_df = (
        pd.concat(snp_dfs).rename({"Start": "Pos"}, axis=1).sort_values(["Chr", "Pos"])
    )
    return snp_df.loc[:, ["Chr", "Pos", "ExonPos", "Ref", "Depth", "Alt", "VAF"]]


# combine SNP data and covData
def combine_CNV(snp_list, cov_list, pon_cov_path="", verbose=False):
    """
    load the coverage_data for a sample and the heteroSNP data and apply the same fullExonCoords
    """

    sample = os.path.basename(snp_list[0].split(".")[0])

    show_output(f"Loading coverage data for sample {sample}")
    cov_df = combine_coverage(
        cov_list,
        pon_cov_path=pon_cov_path,
        verbose=verbose,
    )
    show_output(f"Loading SNP data for sample {sample}")
    snp_df = combine_SNP(snp_list, verbose=verbose)
    show_output(f"Finished loading sample {sample}", color="success")
    return cov_df, snp_df
