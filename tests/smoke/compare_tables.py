#!/usr/bin/env python3
"""Compare two tab-separated tables (plain or .gz) by content.

Same columns and rows (in order); text columns must match exactly; numeric columns within a tolerance.
Usage: compare_tables.py EXPECTED ACTUAL [--atol 1e-6] [--col-atol EB=0.002 ...]
Exit 0 and print PASS, or exit 1 and print what differs (first 10 cells).
"""

import argparse
import sys

import numpy as np
import pandas as pd


def read(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("expected")
    ap.add_argument("actual")
    ap.add_argument("--atol", type=float, default=1e-6)
    ap.add_argument("--col-atol", action="append", default=[], help="COLUMN=TOL")
    a = ap.parse_args()
    col_tol = {k: float(v) for k, v in (x.split("=", 1) for x in a.col_atol)}
    exp, act = read(a.expected), read(a.actual)
    name = a.actual
    if list(exp.columns) != list(act.columns):
        print(f"FAIL {name}: columns differ\n  expected {list(exp.columns)}\n  actual   {list(act.columns)}")
        return 1
    if len(exp) != len(act):
        print(f"FAIL {name}: {len(act)} rows, expected {len(exp)}")
        return 1
    bad = []
    for col in exp.columns:
        e, x = exp[col], act[col]
        en, xn = pd.to_numeric(e, errors="coerce"), pd.to_numeric(x, errors="coerce")
        numeric = en.notna() & xn.notna()
        tol = col_tol.get(col, a.atol)
        diff = (en[numeric] - xn[numeric]).abs()
        for i in diff[diff > tol].index[:10]:
            bad.append(f"  row {i} {col}: expected {e[i]} actual {x[i]} (tol {tol})")
        text = ~numeric
        for i in e[text][(e[text] != x[text])].index[:10]:
            bad.append(f"  row {i} {col}: expected {e[i]!r} actual {x[i]!r}")
    if bad:
        print(f"FAIL {name}: {len(bad)} differing cells (first shown)\n" + "\n".join(bad[:10]))
        return 1
    print(f"PASS {name} ({len(act)} rows x {len(act.columns)} columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
