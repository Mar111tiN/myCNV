import os 
from codeCNV.heteroSNP import get_heteroSNP
from script_utils import set_path, show_output

# get the run_shell function to be passed to running code
# the snakemake object has to be passed to retrieve the proper scripts folder
# run_shell = set_path('codeCNV', snakemake)


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output
    cc = c['CNV']['hetSNP']
    output = s.output.bedCov

    # get the run_shell function to be passed to running code
    # the snakemake object has to be passed to retrieve the proper scripts folder
    run_shell = set_path('codeCNV', s)

    ########## CONFIG #######################
    # squeeze out the config for get_coverage

    config = {
        # paths
        'bedfile': os.path.join(c['paths']['mystatic'], c['ref']['bed_file']),
        'genome_split_path': os.path.join(c['paths']['mystatic'], c['ref']['genome_split']),
        'SNPdb_path': os.path.join(c['paths']['mystatic'], c['ref']['dbsnp_split']),
        # params
        'q': cc['MAPQ'],
        'Q': cc['Q'],
        'minVAF': cc['minVAF'],
        'minDepth': cc['minDepth'],
        'run_shell': run_shell

    }

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
