"""GC and mappability track generators (plan D22): formats the readers accept, values, the public GC track slice."""

import gzip
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from conftest import DATA
from mycnv import trackgen
from mycnv.tracks import load_gc, load_genmap

REF7 = f"{DATA}/ref/chr7_500k.fa.gz"  # GRCh38.GATK chr7:1-500,000 (public reference sequence), bgzipped + .fai
GC7 = f"{DATA}/tracks/gc"
HAS_GENMAP = shutil.which("genmap") is not None
needs_genmap = pytest.mark.skipif(not HAS_GENMAP or shutil.which("samtools") is None,
                                  reason="genmap/samtools not installed (both are in environment.yml)")
pytestmark = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed (in environment.yml)")
MYCNV = [sys.executable, "-m", "mycnv.cli"]


def write_fasta(path, seqs: dict, width=60):
    with open(path, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            fh.write("\n".join(s[i : i + width] for i in range(0, len(s), width)) + "\n")
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    return str(path)


def brute_gc(seq: str):
    up = seq.upper()
    rows = []
    for S in range(0, len(seq) + 10, 10):
        w = up[max(S - 50, 0) : max(min(S + 50, len(seq)), 0)]  # 1-based [S-49, S+50] = 0-based [S-50, S+50)
        gc, at = sum(w.count(b) for b in "GC"), sum(w.count(b) for b in "AT")
        if gc + at:  # numpy rounding (ties at x.xx5 occur only at N-gap edges: 18 on chr21)
            rows.append((S, float(np.round(gc / (gc + at), 2))))
    return pd.DataFrame(rows, columns=["Start", "GC/AT"])


# ---------- GC -------------------------------------------------------------------------------------------------
def test_gc_table_matches_a_direct_window_count():
    rng = np.random.default_rng(1)
    seq = "".join(rng.choice(list("ACGTacgt"), size=1237))
    seq = seq[:300] + "N" * 400 + seq[700:]  # a gap longer than the window: rows inside it are left out
    seq_arr = np.frombuffer(seq.encode(), dtype=np.uint8) & 0xDF
    got = trackgen.gc_table(seq_arr)
    exp = brute_gc(seq)
    assert got["Start"].tolist() == exp["Start"].tolist()
    assert np.allclose(got["GC/AT"], exp["GC/AT"])
    assert not got["Start"].between(360, 640).any()  # windows fully inside the N run have no row


def test_gc_track_reproduces_the_public_track_slice():
    """chr7:1-500,000 of the GATK bundle vs the public static track slice: same rows, 99.96 % equal values (the old
    track was made from a FASTA that differs from the bundle at a few positions)."""
    new = trackgen.gc_table(trackgen.read_sequence(REF7, "chr7"))
    old = pd.read_csv(f"{GC7}/chr7.gc100-10.gz", sep="\t")
    new = new[(new["Start"] >= old["Start"].min()) & (new["Start"] <= old["Start"].max())]
    m = old.merge(new, on="Start", how="outer", suffixes=("_old", "_new"), indicator=True)
    assert (m["_merge"] == "right_only").sum() == 0  # no row the old track lacks
    assert (m["_merge"] == "both").mean() > 0.999
    both = m[m["_merge"] == "both"]
    assert np.isclose(both["GC/AT_old"], both["GC/AT_new"]).mean() > 0.999


def test_make_gc_track_writes_what_load_gc_reads(tmp_path):
    rng = np.random.default_rng(2)
    fa = write_fasta(tmp_path / "g.fa", {"chrA": "".join(rng.choice(list("ACGT"), 900)),
                                         "chrB": "".join(rng.choice(list("ACGT"), 700))})
    out = tmp_path / "gc"
    paths = trackgen.make_gc_track(fa, str(out), threads=2)
    assert sorted(p.split("/")[-1] for p in paths) == ["chrA.gc100-10.gz", "chrB.gc100-10.gz"]
    with gzip.open(out / "chrA.gc100-10.gz", "rt") as fh:
        header, first = fh.readline(), fh.readline()
    assert header == "Start\tGC/AT\n" and len(first.split("\t")[1].strip()) == 4  # 2 decimals, e.g. 0.40
    df = load_gc(str(out), "chrB")
    assert list(df.columns) == ["Pos", "GCratio"] and df["Pos"].iloc[0] == 0


def test_gc_unknown_sequence_and_missing_fai(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", {"chrA": "ACGT" * 50})
    with pytest.raises(ValueError, match="not in"):
        trackgen.read_sequence(fa, "chrZ")
    (tmp_path / "h.fa").write_text(">x\nACGT\n")
    with pytest.raises(FileNotFoundError, match="samtools faidx"):
        trackgen.fasta_chroms(str(tmp_path / "h.fa"))


# ---------- mappability ----------------------------------------------------------------------------------------
def test_parse_modes():
    assert trackgen.parse_modes("30_0,75_1") == [(30, 0), (75, 1)]
    assert trackgen.parse_modes(["100_2"]) == [(100, 2)]
    with pytest.raises(ValueError, match="k-mer length"):
        trackgen.parse_modes("30-0")


def test_change_point_encoding_round_trips():
    v = np.array([1.0, 1.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.333])
    enc = trackgen.encode_change_points(v)
    assert enc.tolist() == [1.0, "", 0.5, "", 0.0, "", 1.0, 0.333]
    dec = pd.Series([np.nan if x == "" else x for x in enc], dtype=float).ffill()
    assert dec.tolist() == v.tolist()
    # make_genmap_track starts every chromosome with a NaN sentinel, so its first row always carries the value
    assert trackgen.encode_change_points(np.array([np.nan, 0.5, 0.5])).tolist() == ["", 0.5, ""]
    # any other NaN would be written as '' and read back as a neighbour's value (review F1, 2026-11-27)
    with pytest.raises(ValueError, match="NaN"):
        trackgen.encode_change_points(np.array([np.nan, 1.0, np.nan, 1.0]))


def test_genmap_positions_without_a_value_are_written_as_0(tmp_path, monkeypatch):
    """genmap leaves no bedgraph interval at N bases and k-mers that reach an N or the contig end. Such positions at
    the start, inside and at the end of the selection read 0 after load_genmap (not a neighbour's value), also across
    chunk borders, and are counted per chromosome and mode."""
    bedgraph = "chrA\t10\t40\t1\nchrA\t40\t60\t0.5\nchrA\t70\t90\t1\n"  # 0-based: no value at 0-9, 60-69, 90-99

    def fake_genmap(*args):
        assert args[0] == "map"
        prefix = args[args.index("-O") + 1]
        with open(prefix + ".bedgraph", "w") as fh:
            fh.write(bedgraph)

    monkeypatch.setattr(trackgen, "_genmap", fake_genmap)
    (tmp_path / "sel.bed").write_text("chrA\t0\t100\n")
    out = tmp_path / "map"
    no_value = {}
    trackgen.make_genmap_track("unused.fa", str(out), chroms=["chrA"], bed=str(tmp_path / "sel.bed"),
                               modes="30_0", index="unused_index", tmp_dir=str(tmp_path), chunk=7,
                               no_value=no_value)
    assert no_value == {"chrA": {"30_0": 30}}
    with gzip.open(out / "hg38_genmap.chrA.txt.gz", "rt") as fh:
        fh.readline()
        assert fh.readline() == "chrA\t1\t0.0\n"  # the first row is written, as 0
    got = load_genmap(str(out), "chrA", modes=["30_0"]).set_index("Pos")["map30_0"]
    exp = [0.0] * 10 + [1.0] * 30 + [0.5] * 20 + [0.0] * 10 + [1.0] * 20 + [0.0] * 10
    assert got.tolist() == exp


def kmer_mappability(seqs: dict, k: int) -> dict:
    """Exact e=0 mappability (both strands, as genmap's default) for every full k-mer start (1-based)."""
    comp = str.maketrans("ACGT", "TGCA")
    counts: dict = {}
    for s in seqs.values():
        for i in range(len(s) - k + 1):
            km = s[i : i + k]
            counts[km] = counts.get(km, 0) + 1
            rc = km.translate(comp)[::-1]
            counts[rc] = counts.get(rc, 0) + 1
    out = {}
    for name, s in seqs.items():
        for i in range(len(s) - k + 1):
            km = s[i : i + k]
            rc = km.translate(comp)[::-1]
            n = counts[km] if km != rc else counts[km] // 2  # a palindrome is one occurrence on both strands
            out[(name, i + 1)] = round(1 / n, 3)
    return out


@pytest.fixture(scope="module")
def repeat_genome(tmp_path_factory):
    d = tmp_path_factory.mktemp("genmap")
    rng = np.random.default_rng(3)
    a = "".join(rng.choice(list("ACGT"), 3000))
    b = "".join(rng.choice(list("ACGT"), 2000))
    rep = a[1000:1120]
    a = a[:2000] + rep + a[2120:]  # the same 120 bp twice on chrA
    comp = str.maketrans("ACGT", "TGCA")
    b = b[:500] + a[200:290].translate(comp)[::-1] + b[590:]  # a reverse-complement copy on chrB
    seqs = {"chrA": a, "chrB": b}
    fa = write_fasta(d / "g.fa", seqs)
    bed = d / "sel.bed"
    bed.write_text("chrA\t150\t400\nchrA\t380\t450\nchrA\t950\t1200\nchrA\t1990\t2150\nchrB\t450\t650\n")
    return d, fa, str(bed), seqs


@needs_genmap
def test_genmap_track_values_layout_and_reader(repeat_genome, tmp_path):
    d, fa, bed, seqs = repeat_genome
    idx = str(tmp_path / "idx")
    trackgen.build_genmap_index(fa, idx)
    out = tmp_path / "map"
    paths = trackgen.make_genmap_track(fa, str(out), bed=bed, modes="30_0,50_0", index=idx, threads=2,
                                       tmp_dir=str(tmp_path))
    assert sorted(p.split("/")[-1] for p in paths) == ["hg38_genmap.chrA.txt.gz", "hg38_genmap.chrB.txt.gz"]
    raw = pd.read_csv(out / "hg38_genmap.chrA.txt.gz", sep="\t")
    assert list(raw.columns) == ["Chr", "Pos", "30_0", "50_0"]
    # dense over the merged selection: chrA 151..450, 951..1200, 1991..2150 (1-based)
    exp_pos = [*range(151, 451), *range(951, 1201), *range(1991, 2151)]
    assert raw["Pos"].tolist() == exp_pos
    assert raw.iloc[0][["30_0", "50_0"]].notna().all()  # the first row carries every value
    assert raw["30_0"].isna().mean() > 0.5  # repeated values are left empty
    dec = load_genmap(str(out), "chrA", modes=["30_0", "50_0"]).set_index("Pos")
    for k, col in ((30, "map30_0"), (50, "map50_0")):
        truth = kmer_mappability(seqs, k)
        pos = [p for p in exp_pos if ("chrA", p) in truth]
        assert dec.loc[pos, col].tolist() == [truth[("chrA", p)] for p in pos], col
    assert (dec.loc[1001:1080, "map30_0"] == 0.5).all()  # inside the duplicated 120 bp
    b = load_genmap(str(out), "chrB", modes=["30_0"]).set_index("Pos")
    assert (b.loc[501:560, "map30_0"] == 0.5).all()  # reverse-complement copy of chrA 201-290


@needs_genmap
def test_genmap_track_n_runs_and_contig_end_read_0(tmp_path):
    """A selected interval that starts and ends inside N runs, and one that reaches the contig end: positions whose
    k-mer holds an N or runs past the end read 0 (e = 0: genmap gives them no value), the rest the exact mappability."""
    rng = np.random.default_rng(4)
    rand = lambda n: "".join(rng.choice(list("ACGT"), n))  # noqa: E731
    s = rand(400) + "N" * 50 + rand(400) + "N" * 30 + rand(200)  # N at 401-450 and 851-880 (1-based), length 1080
    seqs = {"chrN": s}
    fa = write_fasta(tmp_path / "n.fa", seqs)
    bed = tmp_path / "sel.bed"
    bed.write_text("chrN\t420\t860\nchrN\t1000\t1080\n")  # 1-based 421-860 (N at both ends), 1001-1080 (contig end)
    idx = str(tmp_path / "idx")
    trackgen.build_genmap_index(fa, idx)
    out = tmp_path / "map"
    no_value = {}
    trackgen.make_genmap_track(fa, str(out), bed=str(bed), modes="30_0,50_0", index=idx, threads=2,
                               tmp_dir=str(tmp_path), no_value=no_value)
    with gzip.open(out / "hg38_genmap.chrN.txt.gz", "rt") as fh:
        fh.readline()
        assert fh.readline() == "chrN\t421\t0.0\t0.0\n"  # an N base: explicit 0, not back-filled
    dec = load_genmap(str(out), "chrN", modes=["30_0", "50_0"]).set_index("Pos")
    sel = [*range(421, 861), *range(1001, 1081)]
    assert dec.index.tolist() == sel
    for k, col in ((30, "map30_0"), (50, "map50_0")):
        truth = kmer_mappability(seqs, k)
        undefined = [p for p in sel if p + k - 1 > len(s) or "N" in s[p - 1 : p - 1 + k]]
        exp = [0.0 if p in undefined else truth[("chrN", p)] for p in sel]
        assert dec[col].tolist() == exp, col
        assert no_value["chrN"][f"{k}_0"] == len(undefined)
        assert (dec.loc[851:860, col] == 0).all() and (dec.loc[1081 - k + 1 : 1080, col] == 0).all()


@needs_genmap
def test_cli_tracks(repeat_genome, tmp_path):
    d, fa, bed, _ = repeat_genome
    res = subprocess.run([*MYCNV, "tracks", "gc", "--fasta", fa, "-o", tmp_path / "gc", "--chrom", "chrB"],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert [p.name for p in (tmp_path / "gc").iterdir()] == ["chrB.gc100-10.gz"]
    idx = tmp_path / "idx"
    res = subprocess.run([*MYCNV, "tracks", "genmap", "--fasta", fa, "-o", tmp_path / "map", "--bed", bed,
                          "--modes", "30_0", "--index", idx, "--tmpdir", tmp_path, "-t", "2"],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert idx.is_dir()  # --index that did not exist is built and kept
    res = subprocess.run([*MYCNV, "tracks", "genmap", "--fasta", fa, "-o", tmp_path / "map2", "--bed", bed,
                          "--modes", "30_0", "--index", idx, "--chrom", "chrB"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr  # reuses the index
    assert [p.name for p in (tmp_path / "map2").iterdir()] == ["hg38_genmap.chrB.txt.gz"]
    res = subprocess.run([*MYCNV, "tracks", "genmap", "--fasta", fa, "-o", tmp_path / "m3", "--modes", "30x"],
                         capture_output=True, text=True)
    assert res.returncode == 2 and "k-mer length" in res.stderr
