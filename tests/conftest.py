import os

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
LEGACY = os.path.join(DATA, "legacy")
SYNTH = os.path.join(DATA, "synth")
GC = os.path.join(DATA, "tracks", "gc")
GENMAP = os.path.join(DATA, "tracks", "genmap")
# synthetic tracks for synth/multi.bed (chr7, chr21, chrX)
SGC = os.path.join(SYNTH, "tracks", "gc")
SGENMAP = os.path.join(SYNTH, "tracks", "genmap")
STUDY_CONFIG = os.path.join(HERE, "study_config_CNV.yaml")
# Patient-derived fixtures (real PON normals, NHL1 study SNPs) and the legacy outputs made from them. Git-ignored
# (plan D14); recreate with tests/make_fixtures.py + tests/synth.py + tests/legacy_outputs.py from the dev-space data.
RESTRICTED = os.path.join(DATA, "restricted")
RLEGACY = os.path.join(RESTRICTED, "legacy")
needs_restricted = pytest.mark.skipif(
    not os.path.isdir(RLEGACY), reason="restricted (patient-derived) fixtures not present: see tests/README.md"
)


def assert_frames_close(new: pd.DataFrame, old: pd.DataFrame, rtol=1e-9, atol=1e-12, exclude=()):
    """Same columns (order), same length; numeric columns close (NaN == NaN), others equal as strings."""
    assert list(new.columns) == list(old.columns), f"columns differ:\n{list(new.columns)}\n{list(old.columns)}"
    assert len(new) == len(old), f"lengths differ: {len(new)} vs {len(old)}"
    for col in new.columns:
        if col in exclude:
            continue
        a, b = new[col].reset_index(drop=True), old[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            ok = np.isclose(a.to_numpy(float), b.to_numpy(float), rtol=rtol, atol=atol, equal_nan=True)
            assert ok.all(), f"{col}: {int((~ok).sum())} values differ, first at row {int(np.argmin(ok))}: " \
                             f"{a[np.argmin(ok)]} vs {b[np.argmin(ok)]}"
        else:
            assert (a.astype(str) == b.astype(str)).all(), f"{col} differs"


def read(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


@pytest.fixture
def study_config():
    from mycnv.config import load_config

    return load_config(STUDY_CONFIG)


ALT = "[AaCcTtGgDdIi]"  # pile2CNV v2.1.0 alt classes (v1.0.1 lacked the reverse-strand insertion "i")


def expected_het_snps(stream: pd.DataFrame, normal_vaf=(0.2, 0.8), min_depth=15) -> pd.DataFrame:
    """pile2CNV's het-SNP table recomputed from a clean pileup stream (Chr Start Ref Depth1 Read1 Depth2 Read2
    [ExonPos]): normal read has an alt, every sample depth >= min_depth, normal VAF within the range."""
    reads = [c for c in stream.columns if c.startswith("Read")]
    cand = stream.loc[stream["Read1"].str.contains(ALT)].copy()
    keep = np.ones(len(cand), bool)
    out = cand[["Chr", "Start"] + (["ExonPos"] if "ExonPos" in cand.columns else [])].rename(columns={"Start": "Pos"})
    for k, read in enumerate(reads, 1):
        depth = cand[f"Depth{k}"].astype(int)
        alt = cand[read].str.count(ALT)
        keep &= (depth >= min_depth).to_numpy()
        out[f"Depth{k}"] = depth
        out[f"VAF{k}"] = np.where(depth > 0, alt / depth.where(depth > 0, 1), -1)
    keep &= ((out["VAF1"] >= normal_vaf[0]) & (out["VAF1"] <= normal_vaf[1])).to_numpy()
    return out.loc[keep].reset_index(drop=True)


def has_reverse_insertion(stream: pd.DataFrame) -> pd.Series:
    return stream[[c for c in stream.columns if c.startswith("Read")]].apply(lambda s: s.str.contains("i")).any(axis=1)
