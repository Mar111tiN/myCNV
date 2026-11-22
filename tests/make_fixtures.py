"""Cut the small regression fixtures in tests/data/ from the somvar test data (run once; outputs are committed).

Source: the dev-space test data T = genome-bundle-testdata/kumo_pipelines/kumo-variant-calling-nf/somvar/,
synced to a local cache with the same layout (default /home/GOFQR/tools/somvar-dev/data):
  studies/NHL/CNV/NHL1/{raw,data}/NHL1_A-B.*                     (old study outputs, legacy coordinates)
  genome_resources/.../SS_HAEv7/PON/somvar_HAEv7_NovaSeq_v1/      (cov/, snp/; matrix/chr7.pon.gz head)
  static/hg38/{gc100-10/chr7.gc100-10.gz, genmap/split/hg38_genmap.chr7.txt.gz, bed_files/SureSelect/hg38/SS_HAEv7_*}

Usage: python tests/make_fixtures.py [DATA_CACHE]
"""

import gzip
import os
import sys

import pandas as pd

T = sys.argv[1] if len(sys.argv) > 1 else "/home/GOFQR/tools/somvar-dev/data"
PON = os.path.join(T, "genome_resources/Homo_sapiens/GRCh38.GATK/DNAseq/SS_HAEv7/PON/somvar_HAEv7_NovaSeq_v1")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RESTRICTED = os.path.join(OUT, "restricted")  # patient-derived: git-ignored, never committed (plan D14)
MATRIX_ROWS = 4000      # PONmatrix rows of chr7 (positions 35213..~300k)
PON_ROWS = 3000         # PON coverage rows per chromosome
RAW_ROWS = 15000        # NHL1 raw coverage bins per chromosome


def bed_slice(bed: str, chrom: str, end: int, out: str) -> None:
    rows = []
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom and int(f[1]) <= end:
                rows.append("\t".join(f[:3]))
    with open(out, "w") as fh:
        fh.write("\n".join(rows) + "\n")


def main() -> None:
    os.makedirs(os.path.join(RESTRICTED, "pon_mini", "cov"), exist_ok=True)
    os.makedirs(os.path.join(RESTRICTED, "pon_mini", "snp"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "tracks", "gc"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "tracks", "genmap"), exist_ok=True)
    os.makedirs(os.path.join(RESTRICTED, "nhl1"), exist_ok=True)

    # 1. PONmatrix head (chr7) + BED slices covering it
    matrix = os.path.join(PON, "matrix_head", "chr7.pon.head300k.tsv")
    with open(matrix) as fh, gzip.open(os.path.join(RESTRICTED, "chr7.pon.gz"), "wt") as out:
        for i, line in enumerate(fh):
            if i > MATRIX_ROWS:
                break
            out.write(line)
    last = pd.read_csv(os.path.join(RESTRICTED, "chr7.pon.gz"), sep="\t", usecols=["Start"])["Start"].max()
    beds = os.path.join(T, "static/hg38/bed_files/SureSelect/hg38")
    bed_slice(os.path.join(beds, "SS_HAEv7_hg38_Padded.bed"), "chr7", last, os.path.join(OUT, "haev7_padded_chr7.bed"))
    bed_slice(os.path.join(beds, "SS_HAEv7_hg38_Covered.bed"), "chr7", last, os.path.join(OUT, "haev7_covered_chr7.bed"))

    # 2. GC / genmap slices for chr7 over the matrix range (the NHL1 raw tables already carry GCratio/map*)
    raw = pd.read_csv(os.path.join(T, "studies/NHL/CNV/NHL1/raw/NHL1_A-B.cov.gz"), sep="\t")
    raw_parts = [raw.loc[raw["Chr"] == c].head(RAW_ROWS) for c in ["chr7", "chr21"]]
    gc = pd.read_csv(os.path.join(T, "static/hg38/gc100-10/chr7.gc100-10.gz"), sep="\t")
    gc = gc.loc[gc["Start"] <= last + 200]
    gc.to_csv(os.path.join(OUT, "tracks", "gc", "chr7.gc100-10.gz"), sep="\t", index=False, compression="gzip")
    gm = pd.read_csv(os.path.join(T, "static/hg38/genmap/split/hg38_genmap.chr7.txt.gz"), sep="\t")
    gm = gm.loc[gm["Pos"] <= last + 200]
    gm.to_csv(os.path.join(OUT, "tracks", "genmap", "hg38_genmap.chr7.txt.gz"), sep="\t", index=False,
              compression="gzip")

    # 3. mini PON: first rows of cov/ and snp/ for chr7, chr21, chrX
    for chrom in ["chr7", "chr21", "chrX"]:
        cov = pd.read_csv(os.path.join(PON, "cov", f"{chrom}.cov.gz"), sep="\t", nrows=PON_ROWS)
        cov.to_csv(os.path.join(RESTRICTED, "pon_mini", "cov", f"{chrom}.cov.gz"), sep="\t", index=False, compression="gzip")
        end = cov["Pos"].max()
        snp = pd.read_csv(os.path.join(PON, "snp", f"{chrom}.snp.gz"), sep="\t")
        snp = snp.loc[snp["Pos"] <= end]
        snp.to_csv(os.path.join(RESTRICTED, "pon_mini", "snp", f"{chrom}.snp.gz"), sep="\t", index=False, compression="gzip")

    # 4. NHL1 raw slices (legacy coordinates): chr7 and chr21 heads, SNPs in the same ranges
    snp = pd.read_csv(os.path.join(T, "studies/NHL/CNV/NHL1/raw/NHL1_A-B.snp"), sep="\t")
    snp_parts = []
    for part in raw_parts:
        chrom, end = part["Chr"].iloc[0], part["Pos"].max()
        snp_parts.append(snp.loc[(snp["Chr"] == chrom) & (snp["Pos"] <= end)])
    raw_parts[0].to_csv(os.path.join(RESTRICTED, "nhl1", "raw_chr7.cov.gz"), sep="\t", index=False, compression="gzip")
    snp_parts[0].to_csv(os.path.join(RESTRICTED, "nhl1", "raw_chr7.snp"), sep="\t", index=False)
    pd.concat(raw_parts).to_csv(os.path.join(RESTRICTED, "nhl1", "raw_chr7_21.cov.gz"), sep="\t", index=False,
                                compression="gzip")
    pd.concat(snp_parts).to_csv(os.path.join(RESTRICTED, "nhl1", "raw_chr7_21.snp"), sep="\t", index=False)
    print("fixtures written to", OUT)


if __name__ == "__main__":
    main()
