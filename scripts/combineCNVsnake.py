
import os
import sys
import pandas as pd

from script_utils import set_path, show_output
from rollingCNV import rollingCNV
# get the run_shell function to be passed to running code and set PYTHONPATH
run_shell = set_path('codeCNV', snakemake)


def main(s):
    '''
    snakemake wrapper for CNV combine and rolling metrices
    '''

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    cnv_config = c['CNV']['combine']

    sample = w.sample

    # hardcoded
    cnv_path = "cnv"
    cnvPON_path = os.path.join(c['paths']['mystatic'], c['ref']['PONcov'])

    # run the function
    show_output(
        f"Combining CNV data for sample {sample} and computing rolling metrices", time=True)
    cov_df, snp_df, rCov_df, rSNP_df = rollingCNV(sample,
                                                  sample_cnv_path=cnv_path,
                                                  PON_cnv_path=cnvPON_path,
                                                  config=cnv_config)

    show_output(
        f"Rolling metrix computation for sample {sample} finished.", time=True, color='success')
    cov_df.to_csv(o.cov, sep='\t', index=False)
    snp_df.to_csv(o.snp, sep='\t', index=False)
    rCov_df.to_csv(o.rcov, sep='\t', index=False)
    rSNP_df.to_csv(o.rsnp, sep='\t', index=False)
    show_output(
        f"Finished writing output for sample {sample}.", time=True, color='success')


if __name__ == "__main__":
    main(snakemake)