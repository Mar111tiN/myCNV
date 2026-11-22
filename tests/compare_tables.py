"""Compare two tab-separated tables column by column (rows matched on key columns).

  python tests/compare_tables.py <a> <b> [key ...]       default keys: Chr Pos
Prints row counts, rows only in one table and, per shared column, how many values differ (rtol 1e-9) and the
max absolute difference.
"""

import sys

import numpy as np
import pandas as pd

a_path, b_path, *keys = sys.argv[1:]
keys = keys or ["Chr", "Pos"]
a = pd.read_csv(a_path, sep="\t")
b = pd.read_csv(b_path, sep="\t")
print(f"rows: a={len(a)} b={len(b)}; columns equal: {list(a.columns) == list(b.columns)}")
a = a.drop_duplicates(keys)
b = b.drop_duplicates(keys)
m = a.merge(b, on=keys, how="outer", suffixes=("_a", "_b"), indicator=True)
print("rows only in a:", int((m["_merge"] == "left_only").sum()), "| only in b:", int((m["_merge"] == "right_only").sum()))
both = m.loc[m["_merge"] == "both"]
for col in [c for c in a.columns if c in b.columns and c not in keys]:
    x, y = both[f"{col}_a"], both[f"{col}_b"]
    if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
        close = np.isclose(x, y, rtol=1e-9, atol=1e-12, equal_nan=True)
        diff = np.nanmax(np.abs(x.to_numpy(float) - y.to_numpy(float))) if (~close).any() else 0.0
        print(f"{col:22s} differ {int((~close).sum()):8d}  max|d| {diff:.3g}")
    else:
        print(f"{col:22s} differ {int((x.astype(str) != y.astype(str)).sum()):8d}")
