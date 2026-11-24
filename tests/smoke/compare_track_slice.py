#!/usr/bin/env python3
"""Check a rebuilt track against a fixture slice: every Start of the slice must be in the new track with the same value.

Usage: compare_track_slice.py SLICE NEW [--atol 0]   (both tab-separated, .gz ok, first column Start)
Exit 0 and print PASS, or exit 1 and print what differs.
"""

import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slice")
    ap.add_argument("new")
    ap.add_argument("--atol", type=float, default=0.0)
    a = ap.parse_args()
    s = pd.read_csv(a.slice, sep="\t")
    n = pd.read_csv(a.new, sep="\t")
    if list(s.columns) != list(n.columns):
        print(f"FAIL {a.new}: columns {list(n.columns)}, expected {list(s.columns)}")
        return 1
    key, val = s.columns[0], s.columns[1]
    m = s.merge(n, on=key, how="left", suffixes=("_exp", "_new"))
    missing = int(m[f"{val}_new"].isna().sum())
    diff = (m[f"{val}_exp"] - m[f"{val}_new"]).abs()
    bad = int((diff > a.atol).sum())
    if missing or bad:
        print(f"FAIL {a.new}: {missing} slice positions missing, {bad} values differ (> {a.atol})")
        return 1
    print(f"PASS {a.new}: all {len(s)} slice positions present and equal ({len(n)} rows in the new track)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
