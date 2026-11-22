"""Deterministic synthetic inputs (no repo imports, so the legacy env can use it too).

- clean paired pileups (cleanpileup -d format, normal = sample 1) over BED regions
- a synthetic 5-normal PONmatrix over the same positions (for the CLI chain test)
- a synthetic T/N sample in the coordinate system of tests/data/restricted/pon_mini (combine step; restricted,
  because it is derived from real PON members' coverage)
Usage: python tests/synth.py   (writes tests/data/synth/, committed; and restricted/sample_mini/, not committed)
"""

import gzip
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SYNTH = os.path.join(DATA, "synth")
RESTRICTED = os.path.join(DATA, "restricted")


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("browser", "track", "#")) or not line.strip():
                continue
            f = line.split("\t")
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def _reads(rng, depth, alt, alt_letter="A"):
    ref = rng.integers(0, 2, size=depth - alt)
    bases = ["." if r else "," for r in ref]
    for _ in range(alt):
        bases.insert(int(rng.integers(0, len(bases) + 1)), alt_letter if rng.random() < 0.5 else alt_letter.lower())
    return "".join(bases)


def make_pileup(regions, seed=1, flank=3, het_every=17):
    """Clean T/N pileup over BED regions (positions start-flank .. end+flank, 1-based)."""
    rng = np.random.default_rng(seed)
    lines = ["Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2"]
    for chrom, start, end in regions:
        level_t = rng.uniform(0.6, 1.6)
        for pos in range(start - flank, end + flank + 1):
            d1 = int(rng.integers(20, 120))
            d2 = max(1, int(d1 * level_t * rng.uniform(0.8, 1.2)))
            if pos % het_every == 0:
                a1 = int(round(d1 * rng.uniform(0.3, 0.7)))
                a2 = int(round(d2 * rng.uniform(0.1, 0.9)))
            elif pos % 101 == 0:
                a1, a2 = int(d1 * 0.95), int(d2 * 0.95)  # homozygous, filtered by the VAF range
            else:
                a1 = int(rng.random() < 0.02)
                a2 = int(rng.random() < 0.02)
            letter = "ACGTID"[pos % 6]
            lines.append(f"{chrom}\t{pos}\tG\t{d1}\t{_reads(rng, d1, a1, letter)}\t{d2}\t{_reads(rng, d2, a2, letter)}")
    return "\n".join(lines) + "\n"


def make_pon_matrix(regions, n=5, seed=3, flank=3, het_every=23):
    """Synthetic PONmatrix (ebscore pile2count format) for n normals over the same positions as make_pileup."""
    rng = np.random.default_rng(seed)
    lines = ["Chr\tStart\tRef\tA\tG\tC\tT\tI\tD\tDepth"]
    for chrom, start, end in regions:
        for pos in range(start - flank, end + flank + 1):
            plus = rng.integers(10, 60, n)
            minus = rng.integers(10, 60, n)
            cells = {}
            alt = "ACTID"[pos % 5]
            for letter in "AGCTID":
                if letter == alt and pos % het_every == 0:
                    ap, am = (plus * rng.uniform(0.3, 0.6, n)).astype(int), (minus * rng.uniform(0.3, 0.6, n)).astype(int)
                else:
                    ap, am = (rng.random(n) < 0.01).astype(int), np.zeros(n, int)
                cells[letter] = "|".join(map(str, ap)) + "=" + "|".join(map(str, am))
            depth = "|".join(map(str, plus)) + "=" + "|".join(map(str, minus))
            lines.append("\t".join([chrom, str(pos), "G"] + [cells[c] for c in "AGCTID"] + [depth]))
    return "\n".join(lines) + "\n"


def make_sample_from_pon(pon_dir, out_dir, seed=7):
    """T/N coverage and SNP tables (tn2cnv output layout) on the bins/SNPs of a PON (same coordinates)."""
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)
    for chrom in ["chr7", "chr21", "chrX"]:
        cov = pd.read_csv(os.path.join(pon_dir, "cov", f"{chrom}.cov.gz"), sep="\t")
        n = len(cov)
        base = cov[["Cov3", "Cov5"]].to_numpy()
        gain = np.where(np.arange(n) > n // 2, 1.5, 1.0)
        cov1 = np.round(base[:, 0] * rng.uniform(0.9, 1.1, n), 1)
        cov2 = np.round(base[:, 1] * gain * rng.uniform(0.9, 1.1, n), 1)
        out = cov[["Chr", "Pos", "ExonPos", "GCratio", "map30_0", "map50_0", "map75_1", "map100_2"]].copy()
        out["Cov1"], out["Cov2"] = cov1, cov2
        out.to_csv(os.path.join(out_dir, f"mini_T-N.{chrom}.cov.gz"), sep="\t", index=False, compression="gzip")
        snp = pd.read_csv(os.path.join(pon_dir, "snp", f"{chrom}.snp.gz"), sep="\t")
        m = len(snp)
        s = snp[["Chr", "Pos", "ExonPos", "map30_0", "map50_0", "map75_1", "map100_2"]].copy()
        s["Depth1"] = rng.integers(20, 200, m)
        s["VAF1"] = np.round(rng.uniform(0.2, 0.8, m), 6)
        s["Depth2"] = rng.integers(20, 200, m)
        s["VAF2"] = np.round(rng.uniform(0.05, 0.95, m), 6)
        s.to_csv(os.path.join(out_dir, f"mini_T-N.{chrom}.snp"), sep="\t", index=False)


def main():
    os.makedirs(SYNTH, exist_ok=True)
    regions = read_bed(os.path.join(DATA, "haev7_covered_chr7.bed"))[:60]
    with open(os.path.join(SYNTH, "chr7.cpile"), "w") as fh:
        fh.write(make_pileup(regions, seed=1))
    multi = [("chr7", 1000, 1080), ("chr7", 1100, 1104), ("chr7", 1200, 1290), ("chr21", 500, 560),
             ("chr21", 600, 700), ("chrX", 2000, 2033)]
    with open(os.path.join(SYNTH, "multi.bed"), "w") as fh:
        fh.write("browser position chr7:1000-1290\ntrack name=synthetic\n")
        fh.write("".join(f"{c}\t{s}\t{e}\tgene{i}\n" for i, (c, s, e) in enumerate(multi)))
    with open(os.path.join(SYNTH, "multi.cpile"), "w") as fh:
        fh.write(make_pileup(multi, seed=2, het_every=7))
    os.makedirs(os.path.join(SYNTH, "matrix"), exist_ok=True)
    with gzip.open(os.path.join(SYNTH, "matrix", "chr7.pon.gz"), "wt") as fh:
        fh.write(make_pon_matrix(regions))
    # derived from real PON members' coverage: restricted (not committed), only if the restricted fixtures exist
    if os.path.isdir(os.path.join(RESTRICTED, "pon_mini")):
        make_sample_from_pon(os.path.join(RESTRICTED, "pon_mini"), os.path.join(RESTRICTED, "sample_mini"))
    print("synthetic inputs written to", SYNTH)


if __name__ == "__main__":
    main()
