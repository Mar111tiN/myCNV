from io import StringIO
import os
import pandas as pd
from subprocess import PIPE, run
from script_utils import show_output, show_command

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


def get_coverage(bam_file, chrom='', config={}):
    '''
    creates a coverage_df for a bam file on a given chromosome
    '''


    # unwrap the tools with s function
    s = config['run_shell']

    run("pwd", shell=True)
    # the -F 1024 flag is neccessary in order to remove duplicate reads
    view_cmd = f"samtools view -F 1024 -q {config['q']} {bam_file} {chrom}"
    cov_cmd = f"{s('bamCoverage.mawk')} | {s('rollingCoverage.mawk')} {config['rollingWindowSize']} | "
    # the 1 at the end is the option for the filterbed tool to output exonic coords
    cov_cmd += f"{s('filterBed.mawk')} {config['bedfile']} {chrom} 1"
    cmd = f"{view_cmd} | {cov_cmd}"
    show_command(cmd, multi=False)
    cov_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return cov_df
