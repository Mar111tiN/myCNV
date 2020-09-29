import os
from heteroSNP import get_heteroSNP
from script_utils import show_output


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards
    sconfig = s.config
    i = s.input
    output = s.output.bedCov

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    config = {
        # paths
        'bedfile': os.path.join(sconfig['paths']['mystatic'], sconfig['ref']['bed_file']),
        'genome_split_path': os.path.join(sconfig['paths']['mystatic'], sconfig['ref']['genome_split']),
        'SNPdb_path': os.path.join(sconfig['paths']['mystatic'], sconfig['ref']['dbsnp_split']),
        # params
        'q': sconfig['hetero_snp']['MAPQ'],
        'Q': sconfig['hetero_snp']['Q'],
        'minVAF': sconfig['hetero_snp']['minVAF'],
        'minDepth': sconfig['hetero_snp']['minDepth']
    }
    # get the correct path to the tools
    script_path = os.path.join(
        sconfig['snakedir'], sconfig['paths']['scripts'])
    rel_tools = sconfig['scripts']
    # convert relative paths to abs paths using dict comprehension
    tools = {key: os.path.join(
        script_path, rel_tools[key]) for key in rel_tools}
    config.update(tools)

    # input files
    bam_file = i.bam
    show_output(
        f"Detecting SNPs of {w.sample} on chrom {w.chrom}!", time=True)
    # run the coverage tool
    snp_df = get_heteroSNP(bam_file,
                           chrom=w.chrom,
                           config=config,
                           )
    show_output(
        f"Writing heteroSNP of {w.sample} on chrom {w.chrom} to {output}", color="success")
    snp_df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main(snakemake)
