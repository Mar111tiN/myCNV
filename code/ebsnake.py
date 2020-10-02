import os
from ebrun import run_eb


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards
    p = s.params
    static_path = s.config['paths']['mystatic']

    pon_list = os.path.join(static_path, s.config['pon_list'])

    # squeeze out the config for get_coverage
    EBconfig = {
        # paths
        'pon_list': pon_list
    }
    # get the correct path to the tools
    rel_tools = s.config['scripts']
    # convert relative paths to abs paths using dict comprehension
    tools = {key: os.path.join(
        s.config['snakedir'], rel_tools[key]) for key in rel_tools}
    EBconfig.update(tools)

    run_eb(
        table=s.input.table,
        tumor_bam=s.input.tumor_bam,
        output=s.output,
        chrom=w.chrom,
        threads=s.threads,
        EBparams=s.config['EBFilter']['params'],
        EBconfig=EBconfig
    )


if __name__ == "__main__":
    main(snakemake)
