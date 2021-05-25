import os
import pandas as pd
from script_utils_CNV import show_output, cmd2df


def get_rawCNV(
    normal_bam="",
    tumor_bam="",
    tumor_bams=[],
    TN_pileup_file="",  # mpileup from normalbam, tumorbam(s)
    clean_TN_pileup_file="",  # mpileup from normalbam, tumorbam(s) with cleanpileup.mawk already done
    pileup_is_clean=True,
    chrom="",
    config={},
    SNP_output="",
):
    """
    wrapper around CLI chain around the core tool pile2CNV.mawk

    """

    # PARAMS
    # unwrap mawk tools
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    sc = config["hetSNP"]
    cc = config["coverage"]

    # create the basic command and unpack required params
    # ##### FILTERBED
    filter_cmd = f"{mawk('filterBed')} {config['bed_file']} -x -c {chrom}"

    if not SNP_output:
        show_output("Output file for heteroSNP is missing!", color="warning")
        return
    # ##### PILE2CNV
    cnv_cmd = f"{mawk('pile2CNV')} -x -o {SNP_output} -v {sc['normalVAF'][0]} -V {sc['normalVAF'][1]} -d {sc['minDepth']} -c {cc['minCov']}"

    # combine
    cmd = f"{filter_cmd} | {cnv_cmd}"

    if clean_TN_pileup_file:
        if os.path.splitext(clean_TN_pileup_file)[1] == ".gz":
            cmd = f"gunzip < {clean_TN_pileup_file} | {cmd}"
        else:
            cmd = f"cat {clean_TN_pileup_file} | {cmd}"
    else:
        # add the cleanup cmd if
        cmd = f"{mawk('cleanpileup')} -d | {cmd}"
        if TN_pileup_file:
            if os.path.splitext(TN_pileup_file)[1] == ".gz":
                cmd = f"gunzip < {TN_pileup_file} | {cmd}"
            else:
                cmd = f"cat {TN_pileup_file} | {cmd}"
        else:
            # pileup has to be done
            # get the bam files
            if tumor_bams:
                bams = " ".join([normal_bam], tumor_bams)
            else:
                bams = f"{normal_bam} {tumor_bam}"
            # get params from config
            pc = config["pileup"]
            split_genome = os.path.join(config["genome_split_path"], f"{chrom}.fa")
            pileup_cmd = f"samtools mpileup -f {split_genome} -l {config['bed_file']} -r {chrom} -q {pc['MAPQ']} -Q {pc['Q']} {bams}"
            cmd = f"{pileup_cmd} | {cmd}"
    try:
        cov_df = cmd2df(cmd, show=True, multi=False)
        return cov_df
    except:
        show_output("There was an error using shell command", color="warning")
        return cmd


def get_heteroSNP(bam_file, chrom, config):
    """
    LEGACY
    creates a table of all heteroSNP for a given chromosome and the respective VAFs
    only looks at positions from SNP database added as SNPdb
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    split_genome = os.path.join(config["genome_split_path"], f"{chrom}.fa")
    snp_bed = os.path.join(config["SNPdb_path"], f"{config['SNPdb']}.{chrom}.snp")

    pileup_cmd = f"samtools mpileup -f {split_genome} -q {config['MAPQ']} -Q {config['Q']} -r {chrom} -l {snp_bed} {bam_file}"
    snp_cmd = f"{mawk('cleanSNP')} | {mawk('snpVAF')}  {config['minVAF']} | {mawk('filterBed')} {config['bedfile']} -c {chrom} -x"
    cmd = f"{pileup_cmd} | {snp_cmd}"

    snp_df = cmd2df(cmd, show=True, multi=False)

    snp_df = snp_df.loc[
        snp_df["Depth"] > config["minDepth"],
        ["Chr", "Start", "ExonPos", "Ref", "Depth", "Alt", "VAF"],
    ]
    return snp_df


def bam2coverage(bam_file, chrom="", config={}):
    """
    LEGACY
    creates a coverage_df for a bam file on a given chromosome
    wrapper for bam2coverage.mawk

    MAWK INFO
    CLI: *samtools view $bam chr7 |
      + extracts the reads for that chromosome
    *bamCoverage [minCoverage=0] |
      + would be better to have chromCoverage (where is it)
      + would make it more performant
    *rollingCoverage [rollingWindow=100] |

      + every half windowSize a mean coverage is written out
    *filterbed $BED chr7 [writeout exomCoords=1]`

      + filters the output for positions covered by the bedfile
      + filter the output to only exon-spanning rows
    make the command run in memory using stringIO
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
