import os
import re
import pandas as pd

from script_utils import set_path
# get the run_shell function to be passed to running code
run_shell = set_path('codeCNV', snakemake)

from pon_coverage import make_PON_coverage
from script_utils import show_output



def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards

    sconfig = s.config
    pon_config = sconfig['combine_pon_coverage']

    o = s.output

    sample_list = s.input

    # the sample list is that huge blob of files
    # samples should be processed chromosome-wise so chroms and samples have to be filtered out from sample_list
    # easiest to be done using pandas :-)
    sample_df = pd.DataFrame(sample_list, columns=['file'])
    sample_df[['sample', 'Chr']] = sample_df['file'].str.extract(
        r"([^/]+)\.(chr[0-9XY]+)\.")

    sample_list = sample_df['sample'].unique()
    chrom_list = sample_df['Chr'].unique()

    config = {
        'sample_PON_path': s.params.bedCov_path,
        'normCov': pon_config['norm_coverage'],
        'stdFactor': pon_config['std_factor'],
        'verbose_output': False
    }

    # output the file
    show_output(
        f"Combining PON coverage for {len(sample_df.index)} samples", time=True)
    show_output(f"{len(sample_list)} samples detected:")
    full_df, filter_df = make_PON_coverage(
        sample_list, chrom_list=chrom_list, config=config)

    for chrom in chrom_list:
        out_file = os.path.join(s.params.chromCov_path, f"{chrom}")
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
