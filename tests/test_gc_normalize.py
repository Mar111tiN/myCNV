"""GC normalisation (steps D and E): a GC class without autosomal coverage must not produce inf."""

import logging

import numpy as np
import pandas as pd

from mycnv.combine import normalize_gc, normalize_gc_col


def _bins():
    # GC class 0.5: three autosomal bins + one chrX bin. 0.8: one autosomal bin without coverage in Cov1 and one
    # chrX bin with coverage (a class that is rare on the autosomes, e.g. a PON combined on a few chromosomes;
    # chrX is excluded from the class means). 0.3: only on chrX (main drops those rows).
    return pd.DataFrame({
        "Chr": ["chr7", "chr7", "chr21", "chr7", "chrX", "chrX", "chrX"],
        "Pos": [10, 20, 30, 40, 50, 60, 70],
        "GCratio": [0.5, 0.5, 0.5, 0.8, 0.5, 0.3, 0.8],
        "Cov1": [10.0, 20.0, 30.0, 0.0, 5.0, 7.0, 6.0],
        "Cov2": [40.0, 40.0, 40.0, 8.0, 20.0, 9.0, 4.0],
    })


def test_class_scaled_to_mean_100():
    out = normalize_gc_col(_bins(), "Cov2")
    gc5 = out.loc[out["GCratio"] == 0.5]
    assert np.allclose(gc5.loc[gc5["Chr"] != "chrX", "Cov2"], 100)
    assert np.isclose(gc5.loc[gc5["Chr"] == "chrX", "Cov2"].iloc[0], 50)
    gc8 = out.loc[out["GCratio"] == 0.8].set_index("Chr")["Cov2"]
    assert np.isclose(gc8["chr7"], 100) and np.isclose(gc8["chrX"], 50)


def test_zero_mean_class_gives_nan_not_inf(caplog):
    with caplog.at_level(logging.WARNING, logger="mycnv"):
        out = normalize_gc(_bins())
    assert not np.isinf(out[["Cov1", "Cov2"]].to_numpy(float)).any()
    rare = out.loc[out["GCratio"] == 0.8]
    assert len(rare) == 2 and rare["Cov1"].isna().all() and np.allclose(rare["Cov2"], [100, 50])
    # the bins are kept (the bin set must stay the same for every sample); only the chrX-only class is dropped
    assert sorted(out["Pos"]) == [10, 20, 30, 40, 50, 70]
    assert "Cov1" in caplog.text and "mean 0" in caplog.text
