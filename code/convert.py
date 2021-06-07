import os
import re
from script_utils_CNV import show_output


def write_ASCAT(cnv_df, sample="", outpath="", mock_germline=True):
    """
    for a sample XX_A-B it writes into outpath/XX/:
    XX_A-B_A_logr_tumor.tsv
    XX_A-B_B_logr_tumor.tsv
    XX_A-B_A_baf_tumor.tsv
    XX_A-B_B_baf_tumor.tsv
    XX_A-B_baf_normal.tsv
    XX_A-B_logr_normal.tsv


    sample comes in shape sample_tumor-normal
    """
    s = sample.split("_")
    sample_name = s[0]
    tumor_name = s[1].split("-")[0]
    normal_name = s[1].split("-")[1]

    # base_file XX_A-B
    base_file = os.path.join(outpath, sample)

    # change cols
    cnv_df.loc[:, "Chr"] = cnv_df["Chr"].str.replace("chr", "")

    # make mock germline
    mock_df = cnv_df.loc[:, ["Chr", "Pos", "VAF1", "log2ratio1_mean"]].rename(
        dict(
            Chr="chr",
            Pos="pos",
        )
    )
    mock_df.loc[:, "VAF1"] = 0.5
    mock_df.loc[:, "log2ratio1_mean"] = 0

    ##### VAF

    vaf_cols = [col for col in cnv_df.columns if col.startswith("VAF")]

    normal_cols = {"Chr": "chrs", "Pos": "pos", vaf_cols[0]: sample}
    normal_baf_file = f"{base_file}_{normal_name}_baf.tsv"
    cnv_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
        normal_baf_file, sep="\t"
    )

    tumor_cols = {"Chr": "chrs", "Pos": "pos", vaf_cols[1]: sample}
    tumor_baf_file = f"{base_file}_{tumor_name}_baf.tsv"
    cnv_df.loc[:, tumor_cols.keys()].rename(tumor_cols, axis=1).to_csv(
        tumor_baf_file, sep="\t"
    )
    if mock_germline:
        # store the mock_df as fake germline sample
        mock_baf_file = f"{base_file}_germline_baf.tsv"
        mock_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
            mock_baf_file, sep="\t"
        )

    # log2ratio
    log2_pat = re.compile(r"log2ratio[0-9]+_mean$")
    log_cols = [col for col in cnv_df.columns if re.match(log2_pat, col)]

    normal_cols = {"Chr": "chrs", "Pos": "pos", log_cols[0]: sample}
    normal_log_file = f"{base_file}_{normal_name}_logr.tsv"
    cnv_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
        normal_log_file, sep="\t"
    )

    tumor_cols = {"Chr": "chrs", "Pos": "pos", log_cols[1]: sample}
    tumor_log_file = f"{base_file}_{tumor_name}_logr.tsv"
    cnv_df.loc[:, tumor_cols.keys()].rename(tumor_cols, axis=1).to_csv(
        tumor_log_file, sep="\t"
    )
    if mock_germline:
        # store the mock_df as fake germline sample
        mock_log_file = f"{base_file}_germline_logr.tsv"
        mock_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
            mock_log_file, sep="\t"
        )

    show_output(
        f"ASCAT output written to {base_file}_[baf/logr]_[tumor|normal].tsv",
        color="success",
    )
    return cnv_df
