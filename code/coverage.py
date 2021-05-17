import os
import pandas as pd
from script_utils_CNV import cmd2df

# get the coverage on a given chromosome
# filter out duplicates
# CLI: *samtools view $bam chr7 |

#   + extracts the reads for that chromosome
# *bamCoverage [minCoverage=0] |

#   + would be better to have chromCoverage (where is it)
#   + would make it more performant
# *rollingCoverage [rollingWindow=100] |

#   + every half windowSize a mean coverage is written out
# *filterbed $BED chr7 [writeout exomCoords=1]`

#   + filters the output for positions covered by the bedfile
#   + filter the output to only exon-spanning rows
# make the command run in memory using stringIO


def get_coverage(bam_file, chrom="", config={}):
    """
    creates a coverage_df for a bam file on a given chromosome
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    # the -F 1024 flag is neccessary in order to remove duplicate reads
    drop_dups = " -F 1024" if config["drop_duplicates"] else ""
    view_cmd = f"samtools view{drop_dups} -q {config['MAPQ']} {bam_file} {chrom}"
    cov_cmd = f"{mawk('bamCoverage')} | {mawk('rollingCoverage')} {config['rollingWindowSize']} | "
    # the 1 at the end is the option for the filterbed tool to output exonic coords
    cov_cmd += f"{mawk('filterBed')} {config['bedfile']} -c {chrom} -x"
    cmd = f"{view_cmd} | {cov_cmd}"

    cov_df = cmd2df(cmd, show=True, multi=False)
    return cov_df
