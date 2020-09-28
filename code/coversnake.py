import os
from coverage import get_coverage
from script_utils import show_output


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards
    p = s.params
    sconfig = s.config
    output = s.output.bedCov

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    config = {
        'bedfile': os.path.join(sconfig['paths']['mystatic'], sconfig['ref']['bed_file']),
        'rollingWindowSize': sconfig['get_coverage']['rollingWindowSize'],
        'q': sconfig['get_coverage']['MAPQ']
    }

    # get the correct path to the tools
    script_path = os.path.join(
        sconfig['snakedir'], sconfig['paths']['scripts'])
    rel_tools = sconfig['scripts']
    # convert relative paths to abs paths using dict comprehension
    tools = {key: os.path.join(
        script_path, rel_tools[key]) for key in rel_tools}
    config.update(tools)

    bam_file = p.bam
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
