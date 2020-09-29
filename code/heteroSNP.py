from io import StringIO
import os
import pandas as pd
from subprocess import PIPE, run
from script_utils import show_output, show_command


def get_heteroSNP(bam_file, chrom, config):
    '''
    creates a table of all heteroSNP for a given chromosome and the respective VAFs
    '''

    split_genome = os.path.join(config['genome_split_path'], f"{chrom}.fa")
    snp_bed = os.path.join(config['SNPdb_path'], f"db153.{chrom}.snp")
    pileup_cmd = f"samtools mpileup -f {split_genome} -q {config['q']} -Q {config['Q']} -r {chrom} -l {snp_bed} {bam_file}"
    snp_cmd = f"{config['cleanSNP']} | {config['snpVAF']} {config['minVAF']} | {config['filterBed']} {config['bedfile']} {chrom} 1"
    cmd = f"{pileup_cmd} | {snp_cmd}"
    show_command(cmd, multi=False)

    snp_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')

    snp_df = snp_df.loc[snp_df['Depth'] > config['minDepth'], [
        'Chr', '', 'ExonPos', 'Ref', 'Depth', 'Alt', 'VAF']]
    return snp_df
