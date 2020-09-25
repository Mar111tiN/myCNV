from io import StringIO
import os
import pandas as pd
from subprocess import PIPE, run
from script_utils import show_output, show_command

# get the coverage on a given chromosome
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
    # unwrap the tools
    bamCoverage = config['bamCoverage']
    rollingCoverage = config['rollingCoverage']
    filterBed = config['filterBed']

    view_cmd = f"samtools view {bam_file} {chrom}"
    cov_cmd = f"{bamCoverage} | {rollingCoverage} {config['rollingWindowSize']} | "
    # the 1 at the end is the option for the filterbed tool to output exonic coords
    cov_cmd += f"{filterBed} {config['bedfile']} {chrom} 1"
    cmd = f"{view_cmd} | {cov_cmd}"
    show_command(cmd, multi=False)
    cov_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return cov_df
