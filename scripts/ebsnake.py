import os

from script_utils import set_path
# get the run_shell function to be passed to running code and set PYTHONPATH
run_shell = set_path('codeEB', snakemake)

from ebrun import run_eb


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output


    static_path = c['paths']['mystatic']

    pon_list = os.path.join(static_path, c['pon_list'])

    EBconfig = dict(
        pon_list= pon_list,
        run_shell=run_shell # pass the shell_function for the tools
    )

    run_eb(
        table=s.input.table,
        tumor_bam=s.input.tumor_bam,
        output=s.output,
        chrom=w.chrom,
        threads=s.threads,
        EBparams=c['EBFilter']['params'],
        EBconfig=EBconfig
    )


if __name__ == "__main__":
    main(snakemake)
