import os

from script_utils import set_path
# get the run_shell function to be passed to running code and set PYTHONPATH
run_shell = set_path('codeCNV', snakemake)

from coverage import get_coverage
from script_utils import show_output


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output
    cc = c['CNV']['coverage']
    output = s.output.bedCov

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    config = {
        'bedfile': os.path.join(c['paths']['mystatic'], c['ref']['bed_file']),
        'rollingWindowSize': cc['rollingWindowSize'],
        'q': cc['MAPQ'],
        'run_shell': run_shell
    }


    bam_file = i.bam
    show_output(
        f"Calculating coverage of {w.sample} on chrom {w.chrom}!", time=True)
    # run the coverage tool
    cov_df = get_coverage(bam_file,
                          chrom=w.chrom,
                          config=config,
                          )
    show_output(
        f"Writing coverage of {w.sample} on chrom {w.chrom} to {output}", color="success")
    cov_df.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == "__main__":
    main(snakemake)
