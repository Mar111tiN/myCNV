"""Whole-sample rolling step on an old study sample (raw/ -> data/), with the legacy or the new code.

  conda run -n somvar-legacy python tests/study_rolling.py legacy /tmp/myCNV-main <raw.cov.gz> <raw.snp> <outprefix>
  conda run -n somvar        python tests/study_rolling.py new    -              <raw.cov.gz> <raw.snp> <outprefix>
Writes <outprefix>.cnv.snp.gz and <outprefix>.cnv.gz. Compare with tests/compare_tables.py.
"""

import os
import sys

import pandas as pd
import yaml

mode, main_dir, cov_file, snp_file, prefix = sys.argv[1:6]
HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "study_config_CNV.yaml")) as fh:
    cfg = yaml.safe_load(fh)["CNV"]

cov_df = pd.read_csv(cov_file, sep="\t")
snp_df = pd.read_csv(snp_file, sep="\t")
if mode == "legacy":
    sys.path.insert(0, os.path.join(main_dir, "code"))
    from combineCNV import filter_cov
    from rollingCov import rolling_coverage
    from rollingSNP import remergeCNV, rolling_snp

    cov_df = rolling_coverage(filter_cov(cov_df, config=cfg), config=cfg)
    snp_out = rolling_snp(snp_df, cov_df, config=cfg)
    cnv_out = remergeCNV(snp_out, cov_df)
else:
    from mycnv.config import load_config
    from mycnv.pipeline import rolling_step

    snp_out, cnv_out, _ = rolling_step(cov_df, snp_df, load_config(os.path.join(HERE, "study_config_CNV.yaml")))
snp_out.to_csv(f"{prefix}.cnv.snp.gz", sep="\t", index=False, compression="gzip")
cnv_out.to_csv(f"{prefix}.cnv.gz", sep="\t", index=False, compression="gzip")
print("written", prefix)
