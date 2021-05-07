import os
import pandas as pd
from script_utils_CNV import show_output, cmd2df


def get_heteroSNP(bam_file, chrom, config):
    """
    creates a table of all heteroSNP for a given chromosome and the respective VAFs
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    bam_chrom = chrom if config["chrom_with_chr"] else chrom.replace("chr", "")
    split_genome = os.path.join(config["genome_split_path"], f"{chrom}.fa")
    snp_bed = os.path.join(config["SNPdb_path"], f"{config['SNPdb']}.{chrom}.snp")

    pileup_cmd = f"samtools mpileup -f {split_genome} -q {config['q']} -Q {config['Q']} -r {bam_chrom} -l {snp_bed} {bam_file}"
    snp_cmd = f"{mawk('cleanSNP')} | {mawk('snpVAF')}  {config['minVAF']} | {mawk('filterBed')} {config['bedfile']} -c {chrom} -x"
    cmd = f"{pileup_cmd} | {snp_cmd}"

    snp_df = cmd2df(cmd, show=True, multi=False)

    snp_df = snp_df.loc[
        snp_df["Depth"] > config["minDepth"],
        ["Chr", "Start", "ExonPos", "Ref", "Depth", "Alt", "VAF"],
    ]
    return snp_df
