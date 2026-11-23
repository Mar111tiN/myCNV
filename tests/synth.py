"""Deterministic synthetic inputs (no repo imports, so the legacy env can use it too).

- clean paired pileups (cleanpileup -d format, normal = sample 1) over BED regions
- a synthetic 5-normal PONmatrix over the same positions (for the CLI chain test), also for the 3-chromosome BED
- GC and genmap tracks for the 3-chromosome BED (`synth/tracks/`)
- a raw samtools mpileup snippet with read-start/end markers and indels (`synth/raw.mpileup`)
- a synthetic CNV PON (`synth/pon/{cov,snp}/`: 8 normals, 2 of them male; chrX outside the PAR) and a T/N sample in its coordinates
  (`synth/sample/`, tn2cnv layout per chromosome, with a gain, a loss and BAF shifts): the public inputs of the
  PON-combine, combine, rolling and ASCAT parity tests (legacy outputs: tests/legacy_outputs.py)
- a synthetic T/N sample in the coordinate system of tests/data/restricted/pon_mini (combine step; restricted,
  because it is derived from real PON members' coverage)
Usage: python tests/synth.py   (writes tests/data/synth/, committed; and restricted/sample_mini/, not committed)
"""

import gzip
import io
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SYNTH = os.path.join(DATA, "synth")
RESTRICTED = os.path.join(DATA, "restricted")
GZ = {"method": "gzip", "mtime": 0}  # reproducible bytes


def gz_text(path):
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


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


MULTI = [("chr7", 1000, 1080), ("chr7", 1100, 1104), ("chr7", 1200, 1290), ("chr21", 500, 560),
         ("chr21", 600, 700), ("chrX", 2000, 2033)]
SYN_BINS = {"chr7": 1500, "chr21": 700, "chrX": 500}  # synthetic CNV PON: bins per chromosome
# first bin start per chromosome: chrX outside PAR1 (GRCh38 chrX:10,001-2,781,479, two copies in XY), where the
# male normals' halved coverage is realistic and the coverage sex call (mycnv.sex) sees it
SYN_START = {"chr7": 100_000, "chr21": 100_000, "chrX": 3_000_000}
SYN_NORMALS = 8
SYN_MALE = (3, 6)  # 1-based normals with one X
MAP_MODES = ["30_0", "50_0", "75_1", "100_2"]


def make_tracks(regions, out_dir, seed=5):
    """GC (`{chrom}.gc100-10.gz`: Start GC/AT on the 10 bp grid) and genmap (`hg38_genmap.{chrom}.txt.gz`: Chr Pos
    30_0 50_0 75_1 100_2, every base) tracks over 0..max(end)+200 of every BED chromosome."""
    rng = np.random.default_rng(seed)
    os.makedirs(os.path.join(out_dir, "gc"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "genmap"), exist_ok=True)
    ends = {}
    for chrom, _, end in regions:
        ends[chrom] = max(ends.get(chrom, 0), end)
    for chrom, end in ends.items():
        start = np.arange(0, end + 201, 10)
        gc = np.round(np.clip(0.5 + 0.2 * np.sin(start / 170) + rng.normal(0, 0.05, len(start)), 0.1, 0.9), 2)
        pd.DataFrame({"Start": start, "GC/AT": gc}).to_csv(
            os.path.join(out_dir, "gc", f"{chrom}.gc100-10.gz"), sep="\t", index=False, compression=GZ)
        pos = np.arange(1, end + 201)
        gm = pd.DataFrame({"Chr": chrom, "Pos": pos})
        for k, mode in enumerate(MAP_MODES):
            gm[mode] = np.where(rng.random(len(pos)) < 0.1, np.round(rng.uniform(0.05, 0.6, len(pos)), 3), 1.0)
        gm.to_csv(os.path.join(out_dir, "genmap", f"hg38_genmap.{chrom}.txt.gz"), sep="\t", index=False,
                  compression=GZ)


def make_raw_mpileup():
    """Raw `samtools mpileup` text (2 samples) with ^X read starts, $ read ends, +N / -N indels on both strands."""
    rows = [
        ("chr7", 1001, "A", ["6", "^].,.,$.,", "!!!!!!", "5", ".,..,", "!!!!!"]),
        ("chr7", 1002, "C", ["6", ".+2AG,+2ag..,.", "!!!!!!", "5", ".,-1T.,.", "!!!!!"]),
        ("chr7", 1003, "G", ["6", ".,-3tca.*,.", "!!!!!!", "5", "T,t^~T.", "!!!!!"]),
        ("chr21", 501, "T", ["4", "..A$a", "!!!!", "4", ".+1C,+1c,.", "!!!!"]),
    ]
    return "".join("\t".join([c, str(p), r] + f) + "\n" for c, p, r, f in rows)


def _syn_bins(rng, n, start=100_000):
    """Bin starts (Pos) and BED-derived bin ExonPos of n bins in regions of 3..40 bins with gaps between them."""
    pos, exon, region = [], [], []
    p, e, r = start, 0, 0
    while len(pos) < n:
        length = int(rng.integers(3, 41))
        pos += [p + 10 * k for k in range(length)]
        exon += [e + 10 * k + 1 for k in range(length)]
        region += [r] * length
        p += 10 * length + 10 * int(rng.integers(1, 30))
        e += 10 * length
        r += 1
    return np.array(pos[:n]), np.array(exon[:n]), np.array(region[:n])


def make_cnv_pon(out_dir, seed=11):
    """Synthetic CNV PON coverage/SNP tables (pon-coverage layout): cov/{chrom}.cov.gz, snp/{chrom}.snp.gz."""
    rng = np.random.default_rng(seed)
    for sub in ["cov", "snp"]:
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
    size = rng.uniform(0.6, 1.4, SYN_NORMALS)
    for chrom, n in SYN_BINS.items():
        pos, exon, region = _syn_bins(rng, n, start=SYN_START[chrom])
        gc = np.round(np.clip(0.5 + 0.25 * np.sin(np.arange(n) / 60) + rng.normal(0, 0.04, n), 0.08, 0.92), 2)
        capture = rng.uniform(0.5, 1.5, region.max() + 1)[region]
        level = 100 * (0.6 + 0.6 * np.exp(-(((gc - 0.5) / 0.25) ** 2))) * capture
        cov = pd.DataFrame({"Chr": chrom, "Pos": pos, "ExonPos": exon, "GCratio": gc})
        low_map = rng.random(n) < 0.05
        for mode in MAP_MODES:
            cov[f"map{mode}"] = np.where(low_map, np.round(rng.uniform(0.05, 0.6, n), 3), 1.0)
        for i in range(SYN_NORMALS):
            x = 0.5 if (chrom == "chrX" and i + 1 in SYN_MALE) else 1.0
            values = level * size[i] * x * rng.uniform(0.85, 1.15, n)
            outlier = rng.random(n) < 0.01
            values[outlier] *= rng.uniform(2, 4, int(outlier.sum()))
            cov[f"Cov{i + 1}"] = np.round(values, 1)
        cov.to_csv(os.path.join(out_dir, "cov", f"{chrom}.cov.gz"), sep="\t", index=False, compression=GZ)
        # pooled SNPs: about every third bin, a few common (PONVAF > 0.45) ones
        has = rng.random(n) < 0.35
        k = np.flatnonzero(has)
        offset = rng.integers(0, 10, len(k))
        snp = pd.DataFrame({"Chr": chrom, "Pos": pos[k] + offset, "ExonPos": exon[k] + offset})
        for mode in MAP_MODES:
            snp[f"map{mode}"] = cov[f"map{mode}"].to_numpy()[k]
        common = rng.random(len(k)) < 0.05
        snp["VAF"] = np.round(np.where(common, rng.uniform(0.5, 0.95, len(k)), rng.uniform(0.06, 0.4, len(k))), 4)
        snp["Depth"] = rng.integers(100, 2000, len(k))
        snp.to_csv(os.path.join(out_dir, "snp", f"{chrom}.snp.gz"), sep="\t", index=False, compression=GZ)


def make_cnv_sample(pon_dir, out_dir, name="syn_T-N", seed=13):
    """T/N tables in tn2cnv layout ({name}.{chrom}.cov.gz, {name}.{chrom}.snp) on the bins/SNPs of the synthetic PON.

    Tumor: gain (x1.5, BAF 1/3 | 2/3) on chr7 bins 500-899, loss (x0.5, BAF ~0.05 | 0.95) on chr21 bins 200-399.
    SNPs: 80 % of the PON SNPs plus 10 % sample-only ones, a dense cluster of 6 SNPs on chr7 (SNP-density filter),
    and on 3 % of the SNPs map values that differ from the PON's in the 4th decimal (main merged the PON SNPs on the
    map columns too and lost those rows' PONVAF).
    """
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)
    for chrom in SYN_BINS:
        pon = pd.read_csv(os.path.join(pon_dir, "cov", f"{chrom}.cov.gz"), sep="\t")
        n = len(pon)
        level = pon[[f"Cov{i + 1}" for i in range(SYN_NORMALS)]].median(axis=1).to_numpy()
        cn = np.ones(n)
        if chrom == "chr7":
            cn[500:900] = 1.5
        if chrom == "chr21":
            cn[200:400] = 0.5
        out = pon[["Chr", "Pos", "ExonPos", "GCratio"] + [f"map{m}" for m in MAP_MODES]].copy()
        out["Cov1"] = np.round(level * 0.9 * rng.uniform(0.9, 1.1, n), 1)
        out["Cov2"] = np.round(level * 1.2 * cn * rng.uniform(0.9, 1.1, n), 1)
        out.to_csv(os.path.join(out_dir, f"{name}.{chrom}.cov.gz"), sep="\t", index=False, compression=GZ)

        psnp = pd.read_csv(os.path.join(pon_dir, "snp", f"{chrom}.snp.gz"), sep="\t")
        snp = psnp.loc[rng.random(len(psnp)) < 0.8, ["Chr", "Pos", "ExonPos"] + [f"map{m}" for m in MAP_MODES]]
        extra_bins = pon.sample(frac=0.1, random_state=seed)
        extra = pd.DataFrame({"Chr": chrom, "Pos": extra_bins["Pos"].to_numpy() + 5,
                              "ExonPos": extra_bins["ExonPos"].to_numpy() + 5})
        for m in MAP_MODES:
            extra[f"map{m}"] = extra_bins[f"map{m}"].to_numpy()
        parts = [snp, extra]
        if chrom == "chr7":
            b = pon.iloc[100]
            cluster = pd.DataFrame({"Chr": chrom, "Pos": b["Pos"] + np.arange(6), "ExonPos": b["ExonPos"] + np.arange(6)})
            for m in MAP_MODES:
                cluster[f"map{m}"] = 1.0
            parts.append(cluster)
        snp = pd.concat(parts).drop_duplicates(["Chr", "Pos"]).sort_values("Pos").reset_index(drop=True)
        m = len(snp)
        shift = rng.random(m) < 0.03
        snp.loc[shift, "map50_0"] = np.round(snp.loc[shift, "map50_0"] - 0.0001, 4)
        bin_idx = np.searchsorted(pon["Pos"].to_numpy(), snp["Pos"].to_numpy(), side="right") - 1
        snp["Depth1"] = rng.integers(20, 200, m)
        snp["VAF1"] = np.round(np.clip(0.5 + rng.normal(0, 0.06, m), 0.2, 0.8), 6)
        snp["Depth2"] = rng.integers(20, 200, m)
        side = rng.random(m) < 0.5
        base = np.full(m, 0.5)
        base = np.where(cn[bin_idx] > 1, np.where(side, 1 / 3, 2 / 3), base)
        base = np.where(cn[bin_idx] < 1, np.where(side, 0.05, 0.95), base)
        snp["VAF2"] = np.round(np.clip(base + rng.normal(0, 0.04, m), 0.01, 0.99), 6)
        snp.to_csv(os.path.join(out_dir, f"{name}.{chrom}.snp"), sep="\t", index=False)


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
        out.to_csv(os.path.join(out_dir, f"mini_T-N.{chrom}.cov.gz"), sep="\t", index=False, compression=GZ)
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
    multi = MULTI
    with open(os.path.join(SYNTH, "multi.bed"), "w") as fh:
        fh.write("browser position chr7:1000-1290\ntrack name=synthetic\n")
        fh.write("".join(f"{c}\t{s}\t{e}\tgene{i}\n" for i, (c, s, e) in enumerate(multi)))
    with open(os.path.join(SYNTH, "multi.cpile"), "w") as fh:
        fh.write(make_pileup(multi, seed=2, het_every=7))
    os.makedirs(os.path.join(SYNTH, "matrix"), exist_ok=True)
    with gz_text(os.path.join(SYNTH, "matrix", "chr7.pon.gz")) as fh:
        fh.write(make_pon_matrix(regions))
    # 3-chromosome PONmatrix, tracks and raw mpileup for the multi.bed CLI tests
    os.makedirs(os.path.join(SYNTH, "matrix_multi"), exist_ok=True)
    for chrom in ["chr7", "chr21", "chrX"]:
        with gz_text(os.path.join(SYNTH, "matrix_multi", f"{chrom}.pon.gz")) as fh:
            fh.write(make_pon_matrix([r for r in multi if r[0] == chrom], seed=4))
    make_tracks(multi, os.path.join(SYNTH, "tracks"))
    with open(os.path.join(SYNTH, "raw.mpileup"), "w") as fh:
        fh.write(make_raw_mpileup())
    # synthetic CNV PON + sample (public parity inputs of steps D-H)
    make_cnv_pon(os.path.join(SYNTH, "pon"))
    make_cnv_sample(os.path.join(SYNTH, "pon"), os.path.join(SYNTH, "sample"))
    # derived from real PON members' coverage: restricted (not committed), only if the restricted fixtures exist
    if os.path.isdir(os.path.join(RESTRICTED, "pon_mini")):
        make_sample_from_pon(os.path.join(RESTRICTED, "pon_mini"), os.path.join(RESTRICTED, "sample_mini"))
    print("synthetic inputs written to", SYNTH)


if __name__ == "__main__":
    main()
