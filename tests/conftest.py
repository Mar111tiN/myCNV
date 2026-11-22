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
