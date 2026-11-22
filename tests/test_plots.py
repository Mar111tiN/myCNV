"""The optional plots work on tiny tracks (one SNP) and use matplotlib's colormaps API."""

import numpy as np
import pandas as pd

from mycnv.pipeline import cnv_plots, raw_plots


def _tables(n_snp):
    cov = pd.DataFrame({"Chr": ["chr7"] * 30 + ["chr21"] * 30, "Pos": np.r_[np.arange(30), np.arange(30)] * 10,
                        "ExonPos": np.r_[np.arange(30), np.arange(30)] * 10 + 1,
                        "FullExonPos": np.arange(60) * 10 + 1, "log2ratio1": 0.0, "log2ratio2": 0.3,
                        "log2ratio1_mean": 0.0, "log2ratio2_mean": 0.3})
    snp = cov.iloc[:n_snp].assign(VAF1=0.5, VAF2=0.6)
    return cov, snp


def test_plots_with_a_single_snp(tmp_path):
    cov, snp = _tables(1)
    files = cnv_plots(snp, cov.assign(VAF1=np.nan, VAF2=np.nan), "S1_T-N", "T", "N", str(tmp_path))
    files += raw_plots(cov, snp, "S1_T-N", "T", "N", str(tmp_path))
    assert len(files) == 8 and all((tmp_path / f).stat().st_size > 0 for f in files)
