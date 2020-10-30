import yaml
import os
import pandas as pd

# ############ SETUP ##############################
configfile: "config/config_P559CNV.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']

# include helper functions
# include: "includes/io.snk"
include: "includes/utils.snk"

chrom_list = get_chrom_list(config)

# get the samples
sample_sheet = os.path.join(config['snakedir'], config['sample_sheet'])
sample_df = get_bam_files(config['bam_folder'], sample_sheet)
print(sample_df)


# ############ INCLUDES ##############################  
include: "includes/CNV.snk"

# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^/.]+",
    chrom = "(chr)?[0-9XY]+",

# ############## MASTER RULE ##############################################

rule all:
    input:
        # expand("bedCov/{sample}.{chrom}.bedCov", sample=pon_df['sample'], chrom=chrom_list),
        # expand("cnv/{sample}.{chrom}.cov", chrom=chrom_list, sample=sample_df.index),
        # expand("cnv/{sample}.{chrom}.snp", chrom=chrom_list, sample=sample_df.index),
        # expand("cnv/{sample}.{chrom}.snpEB", chrom=chrom_list, sample=sample_df.index),
        expand("CNV/{sample}.cov", sample=sample_df.index),
        expand("CNV/{sample}.snp", sample=sample_df.index),
        expand("CNV/{sample}.roll.cov", sample=sample_df.index),
        expand("CNV/{sample}.roll.snp", sample=sample_df.index)

###########################################################################

# print out of installed tools
onstart:
    print("    PON COVERAGE PIPELINE STARTING.......")
    # write config to the results directory
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")

    # cleanup
    if config['cleanup']:
        print('I would like to cleanup but I don\'t know how!')
