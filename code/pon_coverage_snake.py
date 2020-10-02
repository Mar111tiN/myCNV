import os
from pon_coverage import make_PON_coverage
from script_utils import show_output


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards
    chrom = w.chrom

    sconfig = s.config
    pon_config = sconfig['combine_pon_coverage']

    o = s.output

    sample_list = s.input

    config = {
        'normCov': pon_config['norm_coverage'],
        'minCov': pon_config['min_coverage'],
        'maxMeanSTD': pon_config['max_mean_std'],
        'stdFactor': pon_config['std_factor']
    }

    # output the file
    show_output(f"Combining PON coverage for chrom {chrom}")
    full_df, filter_df = make_PON_coverage(chrom, sample_list, config)

    show_output(
        f"Writing full coverage for chrom {chrom} to {o.fullCov}", color="success")
    full_df.to_csv(o.fullCov, sep='\t', index=False, compression="gzip")
    show_output(
        f"Writing filtered coverage for chrom {chrom} to {o.fullCov}", color="success")
    filter_df.to_csv(o.filterCov, sep='\t', index=False, compression="gzip")


if __name__ == "__main__":
    main(snakemake)
