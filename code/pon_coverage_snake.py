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
        'sample_PON_path': s.params.bedCov_path,
        'normCov': pon_config['norm_coverage'],
        'stdFactor': pon_config['std_factor']
    }

    # output the file
    show_output(f"Combining PON coverage", time=True)
    full_df, filter_df = make_PON_coverage(chrom, sample_list, config)

    for chrom in chrom_list:
        out_file = os.path.join(s.params.chromCov_path, "f{chrom}")
        full_file = f"{out_file}.full.csv.gz"
        show_output(
            f"Writing full coverage for chrom {chrom} to {full_file}", color="success")
        full_df.query('Chr == @chrom').to_csv(
            full_file, sep='\t', index=False, compression="gzip")

        filtered_file = f"{out_file}.filtered.csv.gz"
        show_output(
            f"Writing filtered coverage for chrom {chrom} to {filtered_file}", color="success")
        filter_df.query('Chr == @chrom').to_csv(filtered_file, sep='\t',
                                                index=False, compression="gzip")


if __name__ == "__main__":
    main(snakemake)
