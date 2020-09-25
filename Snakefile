import yaml
import os
import pandas as pd

# ############ SETUP ##############################
configfile: "config/config_devel.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']

# include helper functions
# include: "includes/io.snk"
include: "includes/utils.snk"

chrom_list = get_chrom_list(config)


pon_list = os.path.join(config['paths']['mystatic'], config['pon_list'])

pon_df = get_pon_df(pon_list)
print(pon_df)
# ############ INCLUDES ##############################  
include: "includes/pon_coverage.snk"


# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^/.]+",
    chrom = "(chr)?[0-9XY]+",

# ############## MASTER RULE ##############################################

rule all:
    input:
        expand("bedCov/{sample}.{chrom}.bedCov", sample=pon_df['sample'], chrom=chrom_list)

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
