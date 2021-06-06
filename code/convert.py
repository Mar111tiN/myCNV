import os
import re
import pandas as pd
from script_utils_CNV import show_output


def write_ASCAT(cnv_df, sample="", outpath=""):
    """
    for a sample XX_A-B it writes into outpath/XX/:
    XX_baf_normal.tsv
    XX_baf_tumor.tsv
    XX_logr_normal.tsv
    XX_logr_tumor.tsv

    sample comes in shape sample_tumor-normal
    """
    s = sample.split("_")
    sample_name = s[0]
    tumor_name = s[1].split("-")[0]
    normal_name = s[1].split("-")[1]

    base_file = os.path.join(outpath, f"{sample_name}/{sample_name}")

    cnv_df.loc[:, "Chr"] = cnv_df["Chr"].str.replace("chr", "")
    # VAF
    vaf_cols = [col for col in cnv_df.columns if col.startswith("VAF")]

    normal_cols = {"Chr": "chrs", "Pos": "pos", vaf_cols[0]: sample_name}
    normal_baf_file = f"{base_file}_baf_normal.tsv"
    cnv_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
        normal_baf_file, sep="\t"
    )

    tumor_cols = {"Chr": "chrs", "Pos": "pos", vaf_cols[1]: sample_name}
    tumor_baf_file = f"{base_file}_baf_tumor.tsv"
    cnv_df.loc[:, tumor_cols.keys()].rename(tumor_cols, axis=1).to_csv(
        tumor_baf_file, sep="\t"
    )

    # loglratio
    log2_pat = re.compile(r"log2ratio[0-9]+_mean$")
    log_cols = [col for col in cnv_df.columns if re.match(log2_pat, col)]

    normal_cols = {"Chr": "chrs", "Pos": "pos", log_cols[0]: sample_name}
    normal_log_file = f"{base_file}_logr_normal.tsv"
    cnv_df.loc[:, normal_cols.keys()].rename(normal_cols, axis=1).to_csv(
        normal_log_file, sep="\t"
    )

    tumor_cols = {"Chr": "chrs", "Pos": "pos", log_cols[1]: sample_name}
    tumor_log_file = f"{base_file}_logr_tumor.tsv"
    cnv_df.loc[:, tumor_cols.keys()].rename(tumor_cols, axis=1).to_csv(
        tumor_log_file, sep="\t"
    )
    show_output(
        f"ASCAT output written to {base_file}_[baf/logr]_[tumor|normal].tsv",
        color="success",
    )
    return cnv_df
