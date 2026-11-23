"""Build the GC-ratio and mappability tracks that mycnv reads (`tracks.load_gc`, `tracks.load_genmap`), from a FASTA.

Plan D22: the kit-specific tracks of the studies are lost, so the tracks are regenerated. The same functions serve new
genome-bundle builds (plan L3/L10: `<bundle>/annotation/{gc,genmap}/<name>-<version>/`).

GC track `{chrom}.gc100-10.gz` (columns `Start  GC/AT`):
    For every Start S on a 10 bp grid, GC/AT = (G + C) / (A + C + G + T) over the 100 bp window [S-49, S+50]
    (1-based, both inclusive, clipped to the chromosome; N and other letters are ignored), written with 2 decimals.
    Rows whose window holds no A/C/G/T base are left out. Checked on chr21: this reproduces the old static track for
    every row where the GATK bundle FASTA and the FASTA the old track was made from agree (99.985 %; the old FASTA had
    sequence in ~2 Mb of chr21 that the bundle masks with N).

Mappability track `hg38_genmap.{chrom}.txt.gz` (columns `Chr  Pos  <k>_<e> ...`, e.g. 30_0 50_0 75_1 100_2):
    genmap (`genmap map -K k -E e`, default reverse-complement counting) mappability = 1 / number of occurrences of
    the k-mer starting at Pos with up to e errors, rounded to 3 decimals. Rows are dense over the selection (BED),
    and a cell is left empty when its value equals the value of the previous row (the encoding of the old static
    track; `load_genmap` restores it with forward-fill). The index must be built on the whole genome the reads are
    aligned to (cross-chromosome repeats), which needs ~6 bytes per base of RAM: run a genome-wide build as a kumo
    job or alone on a >= 64 GB host, never next to other work (skill `somvar-port` -> references/genome-tracks.md).
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from .tools import log
from .tracks import GC_PATTERN, GENMAP_MODES, GENMAP_PATTERN

GC_STEP = 10
GC_BLOCKS = 10  # window = 10 blocks of 10 bp: [S-49, S+50] for S = 10 * k
GENMAP_DECIMALS = 3


# ---------------------------------------------------------------- FASTA access ---------------------------------
def fasta_chroms(fasta: str) -> list[str]:
    """Sequence names of an indexed FASTA (plain or bgzipped; needs `<fasta>.fai`, made with `samtools faidx`)."""
    fai = fasta + ".fai"
    if not os.path.isfile(fai):
        raise FileNotFoundError(f"{fai} not found: index the FASTA with `samtools faidx {fasta}`")
    with open(fai) as fh:
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


def read_sequence(fasta: str, chrom: str) -> np.ndarray:
    """One sequence of an indexed FASTA as an upper-case uint8 array (via `samtools faidx`)."""
    if chrom not in fasta_chroms(fasta):
        raise ValueError(f"{chrom} is not in {fasta}")
    out = subprocess.run(["samtools", "faidx", fasta, chrom], check=True, capture_output=True).stdout
    body = out.split(b"\n", 1)[1].replace(b"\n", b"").replace(b"\r", b"")
    return np.frombuffer(body, dtype=np.uint8) & 0xDF  # upper case; N stays N


# ---------------------------------------------------------------- GC -------------------------------------------
def gc_table(seq: np.ndarray) -> pd.DataFrame:
    """GC/AT on the 10 bp grid for one sequence (see the module docstring)."""
    n = len(seq)
    n_blocks = -(-n // GC_STEP)
    pad = np.full(n_blocks * GC_STEP - n, ord("N"), dtype=np.uint8)
    s = np.concatenate([seq, pad]).reshape(n_blocks, GC_STEP)
    gc = ((s == ord("G")) | (s == ord("C"))).sum(axis=1, dtype=np.int32)
    acgt = gc + ((s == ord("A")) | (s == ord("T"))).sum(axis=1, dtype=np.int32)
    cgc = np.concatenate([[0], np.cumsum(gc, dtype=np.int64)])
    cacgt = np.concatenate([[0], np.cumsum(acgt, dtype=np.int64)])
    # Start S = 10 k (k = 0 .. n_blocks): window [10k-49, 10k+50] = blocks k-5 .. k+4 (block b = bases 10b+1 .. 10b+10)
    k = np.arange(n_blocks + 1)
    lo = np.clip(k - GC_BLOCKS // 2, 0, n_blocks)
    hi = np.clip(k + GC_BLOCKS // 2, 0, n_blocks)
    win_gc = cgc[hi] - cgc[lo]
    win_acgt = cacgt[hi] - cacgt[lo]
    keep = win_acgt > 0
    ratio = np.round(win_gc[keep] / win_acgt[keep], 2)
    return pd.DataFrame({"Start": k[keep] * GC_STEP, "GC/AT": ratio})


def _gc_one(args) -> str:
    fasta, chrom, out = args
    df = gc_table(read_sequence(fasta, chrom))
    tmp = out + ".part"
    df.to_csv(tmp, sep="\t", index=False, float_format="%.2f", compression="gzip")
    os.replace(tmp, out)
    log.info("GC track %s: %d rows", os.path.basename(out), len(df))
    return out


def make_gc_track(fasta: str, out_dir: str, chroms=None, threads: int = 1, pattern: str = GC_PATTERN) -> list[str]:
    """Write `<out_dir>/{chrom}.gc100-10.gz` for each chromosome (default: every FASTA sequence). One chromosome
    per process, ~12 bytes per base peak (chr21: 8.6 s, 0.57 GB; chr1: ~3 GB). Returns the written paths."""
    chroms = list(chroms) if chroms else fasta_chroms(fasta)
    os.makedirs(out_dir, exist_ok=True)
    jobs = [(fasta, c, os.path.join(out_dir, pattern.format(chrom=c))) for c in chroms]
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            return list(pool.map(_gc_one, jobs))
    return [_gc_one(j) for j in jobs]


# ---------------------------------------------------------------- mappability ----------------------------------
def parse_modes(modes) -> list[tuple[int, int]]:
    """["30_0", "75_1"] or "30_0,75_1" -> [(30, 0), (75, 1)]"""
    if isinstance(modes, str):
        modes = [m for m in modes.split(",") if m]
    out = []
    for m in modes:
        try:
            k, e = (int(x) for x in str(m).split("_"))
        except ValueError:
            raise ValueError(f"mappability mode {m!r} is not <k-mer length>_<errors>, e.g. 30_0") from None
        out.append((k, e))
    return out


def _genmap(*args) -> None:
    exe = shutil.which("genmap")
    if exe is None:
        raise FileNotFoundError("genmap not found on PATH (bioconda package `genmap`)")
    # --version-check 0: genmap otherwise contacts its update server, which hangs on hosts without internet access
    cmd = [exe, *map(str, args), "--version-check", "0"]
    log.debug("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def build_genmap_index(fasta: str, index_dir: str) -> str:
    """`genmap index` of the whole FASTA into index_dir (must not exist). Needs ~6 bytes per base of RAM."""
    if os.path.exists(index_dir):
        raise FileExistsError(f"{index_dir} exists: pass it as an existing index, or remove it")
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(index_dir))) as tmp:
        src = fasta
        if fasta.endswith((".gz", ".bgz")):  # genmap needs plain text
            src = os.path.join(tmp, "genome.fa")
            with gzip.open(fasta, "rb") as fin, open(src, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        log.info("genmap index of %s -> %s", fasta, index_dir)
        _genmap("index", "-F", src, "-I", index_dir)
    return index_dir


def _selection(fasta: str, chroms: list[str], bed: str | None) -> pd.DataFrame:
    """Selection intervals (0-based, half-open), merged and sorted, restricted to chroms."""
    if bed is None:
        sizes = {}
        with open(fasta + ".fai") as fh:
            for line in fh:
                f = line.split("\t")
                sizes[f[0]] = int(f[1])
        return pd.DataFrame({"chrom": chroms, "start": 0, "end": [sizes[c] for c in chroms]})
    rows = []
    with open(bed) as fh:
        for line in fh:
            if line.startswith(("browser", "track", "#")) or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if f[0] in chroms:
                rows.append((f[0], int(f[1]), int(f[2])))
    sel = pd.DataFrame(rows, columns=["chrom", "start", "end"]).sort_values(["chrom", "start", "end"])
    merged = []
    for chrom, df in sel.groupby("chrom", sort=False):
        cs, ce = None, None
        for s, e in zip(df["start"], df["end"]):
            if cs is not None and s <= ce:
                ce = max(ce, e)
            else:
                if cs is not None:
                    merged.append((chrom, cs, ce))
                cs, ce = s, e
        if cs is not None:
            merged.append((chrom, cs, ce))
    out = pd.DataFrame(merged, columns=["chrom", "start", "end"])
    order = {c: i for i, c in enumerate(chroms)}
    return out.sort_values(["chrom", "start"], key=lambda s: s.map(order) if s.name == "chrom" else s)


def _read_bedgraph(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "value"],
                       dtype={"chrom": str, "start": np.int64, "end": np.int64, "value": np.float64})


def _values_at(bg: pd.DataFrame, positions1: np.ndarray) -> np.ndarray:
    """Value of a (single-chromosome, sorted) bedgraph at 1-based positions; NaN where no interval covers them."""
    starts, ends, vals = bg["start"].to_numpy(), bg["end"].to_numpy(), bg["value"].to_numpy()
    p0 = positions1 - 1
    i = np.searchsorted(starts, p0, side="right") - 1
    ok = (i >= 0) & (p0 < ends[np.clip(i, 0, None)])
    out = np.full(len(p0), np.nan)
    out[ok] = vals[i[ok]]
    return out


def encode_change_points(values: np.ndarray) -> np.ndarray:
    """Object array: the value where it differs from the previous row, '' where it repeats (the first row is filled)."""
    out = np.empty(len(values), dtype=object)
    if not len(values):
        return out
    same = np.zeros(len(values), dtype=bool)
    same[1:] = (values[1:] == values[:-1]) | (np.isnan(values[1:]) & np.isnan(values[:-1]))
    out[:] = ""
    show = ~same & ~np.isnan(values)
    out[np.flatnonzero(show)] = values[show].tolist()  # Python floats: written as 1.0, 0.5, 0.333
    return out


def make_genmap_track(
    fasta: str,
    out_dir: str,
    chroms=None,
    bed: str | None = None,
    modes=GENMAP_MODES,
    index: str | None = None,
    threads: int = 1,
    tmp_dir: str | None = None,
    pattern: str = GENMAP_PATTERN,
    chunk: int = 5_000_000,
) -> list[str]:
    """Write `<out_dir>/hg38_genmap.{chrom}.txt.gz` for each chromosome (default: every FASTA sequence).

    bed: selection (e.g. the kit's padded BED, or exons +-100 bp for a kit-independent exome track); rows are dense
         over it. Without a BED the whole chromosome is written (large: one row per base).
    index: an existing genmap index of the whole genome (`build_genmap_index`); if None, one is built from `fasta`
           into a temp folder and removed afterwards. A genome-wide build needs ~6 bytes per base of RAM.
    threads: genmap -T (genmap's own default is 32).
    """
    chroms = list(chroms) if chroms else fasta_chroms(fasta)
    mode_list = parse_modes(modes)
    names = [f"{k}_{e}" for k, e in mode_list]
    os.makedirs(out_dir, exist_ok=True)
    written = []
    with tempfile.TemporaryDirectory(prefix="mycnv_genmap_", dir=tmp_dir) as tmp:
        idx = index or build_genmap_index(fasta, os.path.join(tmp, "index"))
        sel = _selection(fasta, chroms, bed)
        if sel.empty:
            raise ValueError("the selection holds no region on " + ", ".join(chroms))
        sel_bed = os.path.join(tmp, "selection.bed")
        sel.to_csv(sel_bed, sep="\t", header=False, index=False)
        graphs = {}
        for (k, e), name in zip(mode_list, names):
            prefix = os.path.join(tmp, f"map_{name}")
            log.info("genmap map -K %d -E %d (%d threads)", k, e, threads)
            _genmap("map", "-K", k, "-E", e, "-I", idx, "-O", prefix, "-S", sel_bed, "-bg", "-T", threads)
            graphs[name] = _read_bedgraph(_find_output(prefix, tmp))
        for chrom in chroms:
            regions = sel[sel["chrom"] == chrom]
            if regions.empty:
                log.warning("no selected region on %s: no mappability track written", chrom)
                continue
            out = os.path.join(out_dir, pattern.format(chrom=chrom))
            part = out + ".part"
            bgs = {n: g[g["chrom"] == chrom].sort_values("start") for n, g in graphs.items()}
            prev = {n: np.array([np.nan]) for n in names}  # NaN before the first row: it carries every value
            with gzip.open(part, "wt") as fh:
                fh.write("\t".join(["Chr", "Pos", *names]) + "\n")
                for pos in _position_chunks(regions, chunk):
                    cols = []
                    for n in names:
                        v = np.round(_values_at(bgs[n], pos), GENMAP_DECIMALS)
                        cols.append(encode_change_points(np.concatenate([prev[n], v]))[1:])
                        prev[n] = v[-1:]
                    df = pd.DataFrame({"Chr": chrom, "Pos": pos, **dict(zip(names, cols))})
                    df.to_csv(fh, sep="\t", header=False, index=False)
            os.replace(part, out)
            written.append(out)
            log.info("mappability track %s written", os.path.basename(out))
    return written


def _find_output(prefix: str, folder: str) -> str:
    """genmap names its bedgraph after -O (single FASTA) or after the FASTA inside -O (several)."""
    for cand in (prefix + ".bedgraph", prefix):
        if os.path.isfile(cand):
            return cand
    hits = [os.path.join(folder, f) for f in os.listdir(folder) if f.startswith(os.path.basename(prefix))
            and f.endswith(".bedgraph")]
    if os.path.isdir(prefix):
        hits += [os.path.join(prefix, f) for f in os.listdir(prefix) if f.endswith(".bedgraph")]
    if len(hits) != 1:
        raise FileNotFoundError(f"genmap bedgraph for {prefix} not found ({hits})")
    return hits[0]


def _position_chunks(regions: pd.DataFrame, chunk: int):
    """1-based positions of the regions ([start, end) 0-based -> start+1 .. end), in chunks of <= chunk."""
    buf, size = [], 0
    for s, e in zip(regions["start"], regions["end"]):
        s0 = s
        while s0 < e:
            take = min(e - s0, chunk - size)
            buf.append(np.arange(s0 + 1, s0 + take + 1, dtype=np.int64))
            size += take
            s0 += take
            if size >= chunk:
                yield np.concatenate(buf)
                buf, size = [], 0
    if buf:
        yield np.concatenate(buf)
