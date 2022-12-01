import os
import pandas as pd
from script_utils_CNV import show_output, cmd2df


def addGCratio(cov_df, chrom="", gc_path="", mode="100-10"):
    """
    adds the GCdata to the coverage data
    gc_path is path to chrom-split gc-data
    """

    # load genmap data for that chromosome
    gc_file = os.path.join(gc_path, f"{chrom}.gc{mode}.gz")
    if os.path.isfile(gc_file):
        try:
            show_output(f"Loading GC data for {chrom} from {gc_file}")
            gc_df = pd.read_csv(gc_file, sep="\t", compression="gzip").rename(
                {"GC/AT": "GCratio", "Start": "Pos"}, axis=1
            )
        except Exception as e:
            show_output(
                f"{e}: Could not load GC data for {chrom} from {gc_file}",
                color="warning",
            )
            return
        # gc data has columns Start -> GC/AT
        if "Chr" in gc_df.columns:
            cov_df = cov_df.merge(gc_df, on=["Chr", "Pos"], how="left")
        else:
            cov_df = cov_df.merge(gc_df, on="Pos", how="left")
        return cov_df
    else:
        show_output(f"Could not find GC file {gc_file}", color="warning")


def addGenmap(*dfs, chrom="", genmap_path="", modes=["30_0", "50_0", "75_1", "100_2"]):
    """
    adds the genmap data to the coverage data
    selects only the columns that are given in modes
    choose from [
        '30_0', '30_1', '30_2',
        '50_0', '50_1', '50_2',
        '75_0','75_1', '75_2', '75_3',
        '100_0', '100_1', '100_2', '100_4',
        '150_0','150_1', '150_2', '150_4']

    """

    # load genmap data for that chromosome
    genmap_file = os.path.join(genmap_path, f"hg38_genmap.{chrom}.txt.gz")

    if os.path.isfile(genmap_file):
        try:
            cols = ["Chr", "Pos"] + modes
            genmap_df = (
                pd.read_csv(genmap_file, sep="\t", compression="gzip")
                .loc[:, cols]
                .fillna(method="ffill").fillna(method="bfill")
            )
            # rename the mappability cols
            mode_rename = {mode: "map" + mode for mode in modes}
            genmap_df = genmap_df.rename(mode_rename, axis=1)
            show_output(f"Loading mappability data for {chrom} from {genmap_file}")
        except Exception as e:
            show_output(
                f"{e}: Could not load mappability data for {chrom} from {genmap_file}",
                color="warning",
            )
            return
        map_dfs = []
        for df in dfs:
            df = df.merge(genmap_df, on=["Chr", "Pos"], how="left")
            # save cols
            cols = df.columns
            base_cols = [
                col
                for col in cols
                if col in ["Chr", "Start", "Pos", "ExonPos", "GCratio"]
                or col.startswith("map")
            ]
            data_cols = [col for col in cols if col not in base_cols]
            df = df.loc[:, base_cols + data_cols]
            map_dfs.append(df)
        return map_dfs
    else:
        show_output(f"Could not find genmap file {genmap_file}", color="warning")


def TN2CNV(
    normal_bam="",
    tumor_bam="",
    tumor_bams=[],
    TN_pileup_file="",  # mpileup from normalbam, tumorbam(s)
    clean_TN_pileup_file="",  # mpileup from normalbam, tumorbam(s) with cleanpileup.mawk already done
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
    except Exception as e:
        show_output(f"There was an error using shell command <<{e}>>", color="warning")
        return cmd

    # add GC
    if "gc_split_path" in config and os.path.isdir(gc_path := config["gc_split_path"]):
        cov_df = addGCratio(cov_df, chrom=chrom, gc_path=gc_path)
    else:
        show_output(f"Could not find GC path {gc_path}", color="warning")

    # add genmap data to both cov and snp data
    if "genmap_split_path" in config and os.path.isdir(
        genmap_path := config["genmap_split_path"]
    ):
        # reload snp_df from temp file
        show_output(f"Reloading raw heteroSNP data from {SNP_output}")
        snp_df = pd.read_csv(SNP_output, sep="\t")
        cov_df, snp_df = addGenmap(cov_df, snp_df, chrom=chrom, genmap_path=genmap_path)
        # resave snp_df
        show_output(f"Resaving annotated heteroSNP data to {SNP_output}")
        snp_df.to_csv(SNP_output, index=False, sep="\t")
    else:
        show_output(f"Could not find genmap path {genmap_path}", color="warning")
    return cov_df, snp_df


def PON2CNV(chrom="", config={}):
    """
    wrapper around CLI chain around the core tool PON2CNV.mawk

    """

    # PARAMS
    # unwrap mawk tools
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    pon_path = config["PON_path"]
    c = config["PONcoverage"]

    # ####BUILD COMMAND #########
    # ### READ PONMATRIX
    matrix_file = os.path.join(pon_path, f"matrix/{chrom}.pon.gz")
    if not os.path.isfile(matrix_file):
        show_output(
            f"PON matrix file {matrix_file} not found! Exiting.", color="warning"
        )
        return
    read_cmd = f"gunzip < {matrix_file}"

    # ### FILTERBED
    filter_cmd = f"{mawk('filterBed')} {config['bed_file']} -x -c {chrom}"

    # ##### PON2CNV
    SNP_file = os.path.join(pon_path, f"snp/{chrom}.snp")

    cnv_cmd = f"{mawk('PON2CNV')} -x -o {SNP_file} -v {c['minVAF']} -d {c['minDepth']} -c {c['minCov']}"

    # combine
    cmd = f"{read_cmd} | {filter_cmd} | {cnv_cmd}"

    try:
        cov_df = cmd2df(cmd, show=True, multi=False)
    except Exception as e:
        show_output(f"There was an error using shell command <<{e}>>", color="warning")
        return cmd

    # add GC
    if "gc_split_path" in config and os.path.isdir(gc_path := config["gc_split_path"]):
        cov_df = addGCratio(cov_df, chrom=chrom, gc_path=gc_path)
    else:
        show_output(f"Could not find GC path {gc_path}", color="warning")

    # add genmap data to both cov and snp data
    if "genmap_split_path" in config and os.path.isdir(
        genmap_path := config["genmap_split_path"]
    ):
        # reload snp_df from temp file
        show_output(f"Reloading PONSNP data from {SNP_file}")
        snp_df = pd.read_csv(SNP_file, sep="\t")
        cov_df, snp_df = addGenmap(cov_df, snp_df, chrom=chrom, genmap_path=genmap_path)
        # resave snp_df
        show_output(f"Resaving annotated heteroSNP data to {SNP_file}.gz")
        snp_df.to_csv(f"{SNP_file}.gz", index=False, sep="\t", compression="gzip")
        os.remove(SNP_file)
    else:
        show_output(f"Could not find genmap path {genmap_path}", color="warning")
    return cov_df, snp_df


# ############## LEGACY


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
