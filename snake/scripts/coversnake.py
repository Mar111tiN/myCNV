import os

from script_utils import show_output, set_path
from codeCNV.coverage import get_coverage


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    # get the run_shell function to be passed to running code and set PYTHONPATH
    run_shell = set_path('codeCNV', s)

    cc = c['CNV']['coverage']
    output = o.bedCov
    # get the run_shell function to be passed to running code
    # the snakemake object has to be passed to retrieve the proper scripts folder
    run_shell = set_path('codeCNV', s)

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    config = {
        'bedfile': os.path.join(c['paths']['mystatic'], c['ref']['bed_file']),
        'rollingWindowSize': cc['rollingWindowSize'],
        'q': cc['MAPQ'],
        'run_shell': run_shell,
        'chrom_with_chr': c['chrom_with_chr'],
        'drop_duplicates': cc['drop_duplicates']
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
