import os
import pandas as pd
from os import system as shell
from ebutils import get_pon_bases, get_sample_pos, compute_matrix2EB_multi, compute_AB2EB_multi
from script_utils import show_output, show_command, run_cmd


def run_eb(table, tumor_bam, output, chrom, threads, EBparams, EBconfig):
    '''
    master function to start eb_computation
    '''

    # s is helper function returning absolute paths to shell tools
    s = EBconfig['run_shell']

    # unpack the paths
    pon_list = EBconfig['pon_list']

    # ############## LOAD DATA ###############################
    show_output(f"Computing EBscore for chrom {chrom}", color='normal')

    # convert the snp file to a df for EB computation:

    # get the sceleton mutation file for that chromosome

    # convert if it is a snp file
    if "snp" in table:
        mut_df = pd.read_csv(table, sep='\t')
        mut_df['Alt'] = mut_df['Alt'].str.extract(r"([ACTGIiDdacgt-])")
        mut_df['End'] = mut_df['Start']
        # convert Insertions
        mut_df.loc[mut_df['Alt'] == "I", 'Ref'] = '-'
        mut_df.loc[mut_df['Alt'] == "I",
                   'Start'] = mut_df.loc[mut_df['Alt'] == "I", 'Start'] + 1
        # convert Deletions
        mut_df.loc[mut_df['Alt'] == "D", 'Alt'] = '-'
        mut_df.loc[mut_df['Alt'] == "D",
                   'Start'] = mut_df.loc[mut_df['Alt'] == "D", 'Start'] - 1

        mut_df = mut_df.loc[:, ['Chr', 'Start', 'End', 'Ref', 'Alt']]
        # set base_name for intermediate files
        base_file = output[0].replace(".snpEB", "")
    # standard file from varscan pipeline
    else:
        mut_df = pd.read_csv(table, sep='\t', index_col=False, header=None, names=[
            'Chr', 'Start', 'End', 'Ref', 'Alt', 'somatic_status', 'TR1', 'TR1+', 'TR2', 'TR2+', 'NR1', 'NR1+', 'NR2', 'NR2+', 'somaticP', 'variantP']).query('Chr == @chrom').iloc[:, :5]
        # set base_name for intermediate files
        base_file = output[0].replace(".EB", "")

    mut_cols = list(mut_df.columns)

    # ############## PILEUP --> MATRIX FILE ##################

    matrix_file = f"{base_file}.matrix"

    # short_cut for debugging
    if os.path.isfile(matrix_file):
        show_output(
            f"Found pileup up {chrom} of {tumor_bam}. Using it..", color='normal')
    else:
        # bed file can contain all chromosomes because chrom restriction comes with the -r parameter
        bed_file = f"{base_file}.bed"
        # create the bed file for mpileup
        shell(f"{s('csv2bed.mawk')} < {table} > {bed_file}")

        # create the pon_list containing the tumor-bam as first file
        sample_list = f"{base_file}.pon"
        # makeponlist removes the sample itself from list if it is part of PoN
        shell(f"{s('makeponlist.sh')} {tumor_bam} {pon_list} {sample_list}")

        show_output(
            f"Piling up {chrom} of {tumor_bam} with Pon List.", color='normal')
        shell(f"cat {sample_list}")
        # do the pileup into the matrix file

        pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']}"
        pileup_cmd += f" -l {bed_file} -r {chrom} -b {sample_list}"
        # cut -f $({pon2cols}< {sample_list}) creates a cut command only including the desired

        pipe_cmd = f"{pileup_cmd} | cut -f $({s('pon2cols.mawk')} < {sample_list}) | {s('cleanpileup.mawk')}| {s('pile2count.mawk')} > {matrix_file}"
        # do the pileup to matrix_file
        show_command(pipe_cmd, multi=False)
        shell(pipe_cmd)
        # cleanup
        shell(f"rm {bed_file} {sample_list}")

        # check if matrix_file has input
        if not os.path.getsize(matrix_file):
            # create empty file
            open(output[0], 'a').close()
            show_output(
                f"Pileup for {chrom} of {tumor_bam} was empty! Created empty file {output[0]}", color='warning')
            return
        # if no error:
        show_output(
            f"Pileup matrix for chrom {chrom} of {tumor_bam} completed.", color='normal')
    # ################ MERGE INTO MUTFILE ######################
    # change mutation positions for deletions in mutation file
    mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] - 1
    # read matrix file into df
    matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=False)
    # merge
    mut_matrix = mut_df.merge(matrix_df, on=['Chr', 'Start'], how='inner')
    # reset deletion positions
    mut_matrix.loc[mut_matrix['Alt'] == "-",
                   'Start'] = mut_matrix['Start'] + 1

    # ####### if using matrix2EBinput.mawk #######################
    # write to file
    mutmatrix_file = f"{base_file}.mutmatrix"
    mut_matrix.to_csv(mutmatrix_file, sep='\t', index=False)

    # convert mutmatrix to direct EBinput
    EB_matrix_input_file = f"{base_file}.EB.matrix"
    shell(f"cat {mutmatrix_file} | {s('matrix2EBinput.mawk')} > {EB_matrix_input_file}")

    # load in the EB.matrix file
    eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')

    # multithreaded computation
    EB_df = compute_matrix2EB_multi(
        eb_matrix, EBparams['fitting_penalty'], threads)

    # add EBscore to columns
    mut_cols.append('EBscore')

    # get the pon_matrix containing the Pon coverages in Alt and Ref
    pon_matrix = get_pon_bases(eb_matrix)
    # transfer PoN-Ref and PoN-Alt to EB_df
    EB_df[['PoN-Ref', 'PoN-Alt']] = pon_matrix[['PoN-Ref', 'PoN-Alt']]
    mut_cols += ['PoN-Ref', 'PoN-Alt']

    # rm unnecessary columns
    EB_df = EB_df[mut_cols]

    # ######### WRITE TO FILE ##############################################

    EB_file = output[0]
    EB_df.to_csv(EB_file, sep='\t', index=False)

    # cleanup
    shell(f"rm {matrix_file} {EB_matrix_input_file}")  # {mutmatrix_file}
    show_output(
        f"Created EBscore for chrom {chrom} of {tumor_bam} and written to {output[0]}", color='success')
