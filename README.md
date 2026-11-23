# mycnv

PON-based copy-number profiles from a paired (normal, tumor) pileup: 10 bp coverage bins as log2ratio vs a Panel of
Normals, het-SNP B-allele frequencies, GC/mappability filters, rolling windows, and ASCAT input. This is the CNV core of
the somVarWES pipeline, packaged for the kumo-variant-calling-nf port (branch `kumodev`; `main` is the legacy code).

## Install

```bash
conda env create -f environment.yml -n mycnv        # unpinned conda-forge/bioconda env (mawk, pandas, matplotlib, ...)
conda run -n mycnv python -m pip install --no-deps --no-build-isolation .
conda run -n mycnv python -m pytest -q              # public tests incl. parity with the legacy code
```
Container: `docker build -t mycnv:<ver> .` (named stages `build` → `test` → `runtime`; `--target test` runs the public
tests). Only `pyproject.toml`, `README.md`, `src/` and `tests/` enter the image; `tests/data/restricted/` never does.
CI (`.github/workflows/ci.yml`, plan D24): on push/PR to `kumodev`, pytest, the Docker test stage, the runtime image
with a contract check, and Trivy 0.69.3: a Dockerfile config scan, the vulnerability gate (library packages, fixed
vulnerabilities, HIGH/CRITICAL; accepted findings go into `.trivyignore` with a reason) and an informational report
that includes the OS packages.

## Inputs

- **Paired pileup** (normal = sample 1): `samtools mpileup -f ref.fa -q 20 -Q 25 -l <padded BED> NORMAL TUMOR |
  cleanpileup -d` (as written by `ebscore pileup`), plain or gzipped. `--raw-pileup` accepts raw mpileup text.
- **PONmatrix** folder (`--ponmatrix`): `matrix/{chrom}.pon.gz` (`Chr Start Ref A G C T I D Depth`, cells
  `s1|…|sN=s1|…|sN`) from `ebscore pon-matrix`, shared with ebscore. mycnv adds `cov/`, `snp/` and `CNV/`, each in a
  folder you name with `-o/--outdir` (it never writes into an input folder).
- **BED** (`--bed`): the kit's **padded** BED, the same file for the PON and the samples (plan D20; overlapping regions are merged, never double-counted). Every chromosome
  must be one block sorted by start (`sort -k1,1 -k2,2n`); other BEDs are rejected.
- **GC** (`--gc`, required): directory with `{chrom}.gc100-10.gz` (`Start GC/AT`), one file with a `Chr` column, or
  one chromosome's file (its name must name the chromosome).
- **Mappability** (`--genmap`, required for samples): directory with `hg38_genmap.{chrom}.txt.gz`
  (`Chr Pos 30_0 …`) or one file. The generic `genmap/split/` track is not the kit-specific `HAEv7split` the studies
  used (see CHANGELOG, "Known differences"). Both tracks can be rebuilt with `mycnv tracks` (next section).

## Tracks (GC and mappability, for genome-bundle builds)

```bash
mycnv tracks gc     --fasta genome.fa -o gc100-10/ [--chrom chr21 ...] [-t 8]
mycnv tracks genmap --fasta genome.fa -o genmap/ --bed kit.padded.bed [--modes 30_0,50_0,75_1,100_2] \
    [--index genmap_index/] [--chrom chr21 ...] [-t 8] [--tmpdir scratch/]
```
- **GC** (`{chrom}.gc100-10.gz`, `Start GC/AT`): GC/(GC+AT) over [S−49, S+50] on a 10 bp grid, 2 decimals, rows
  without any A/C/G/T base left out. It reproduces the old static track wherever the FASTAs agree (chr21: 99.985 % of
  rows). ~12 bytes per base peak per process (chr21 8.6 s / 0.57 GB), so it is safe genome-wide.
- **Mappability** (`hg38_genmap.{chrom}.txt.gz`): genmap mappability (1 / occurrences on both strands, ≤ e errors)
  of the k-mer starting at each position, 3 decimals, dense over `--bed`, a cell left empty when it repeats the
  previous row (the old encoding). The index must cover the **whole genome** the reads are aligned to (cross-chromosome
  repeats) and needs ~6 bytes per base of RAM: build it once with `--index DIR` (kept and reused). A whole-genome
  index is a job for a kumo job or a host where nothing else runs with ≥ 64 GB; never next to other work.
- Python: `mycnv.trackgen.make_gc_track`, `make_genmap_track`, `build_genmap_index`, `gc_table`.

## Commands

```bash
# CNV PON (once per kit + normal set): coverage/SNPs per chromosome, then the combined PON tables
mycnv pon-coverage --ponmatrix PONmatrix --bed kit.padded.bed --gc gc100-10/ [--genmap genmap/split/] \
    -o ponmatrix_cnv -t 8                               # -> ponmatrix_cnv/cov/{chrom}.cov.{gz,json}, snp/{chrom}.snp.gz
mycnv pon-combine --ponmatrix ponmatrix_cnv -o ponmatrix_cnv   # -> CNV/pon.cov.{full,fullGC,filter,filterGC}.gz,
                                                               #    CNV/pon.snp.gz, CNV/pon.json
# per chromosome (scatter): pon-coverage --matrix PONmatrix/matrix/chr7.pon.gz ... -o cov_chr7  (repeat --chrom/--matrix)

# per pair, all steps in one task (always the whole sample)
mycnv sample --pileup pair.pileup.gz --bed kit.padded.bed --gc gc100-10/ --genmap genmap/split/ --ponmatrix ponmatrix_cnv \
    --sample P1_T-N -o CNV/P1 [--tumor-name T --normal-name N] [--plots]
# -> CNV/P1/raw/P1_T-N.{cov.gz,snp}  data/P1_T-N.{cnv.snp.gz,cnv.gz}  pre/P1_T-N_{N,T,germline}_{baf,logr}.tsv  [plots/]

# or step by step (same outputs); a per-chromosome split stops after tn2cnv, `combine` takes all its outputs
mycnv tn2cnv   --pileup pair.pileup.gz --bed kit.padded.bed --gc … --genmap … --out-cov s.chr7.cov.gz --out-snp s.chr7.snp \
    [--chrom chr7 ...]
mycnv combine  --cov s.chr*.cov.gz --snp s.chr*.snp --ponmatrix ponmatrix_cnv --out-cov raw.cov.gz --out-snp raw.snp
mycnv rolling  --cov raw.cov.gz --snp raw.snp --out-snp data.cnv.snp.gz --out-cnv data.cnv.gz
mycnv ascat-input --cnv data.cnv.snp.gz --sample P1_T-N -o pre/
mycnv config   # print the effective parameters (defaults = the somVarWES study config)
```
Every command takes `--config my.yaml` (structure of `mycnv config`; unknown keys are reported and ignored),
`-t/--threads` (parallel only in `pon-coverage`), `-v/--verbose` (prints the shell pipelines) and `--quiet`. Commands
that run the mawk tools take `--tmpdir`. `--chrom` is repeatable; `--cov`, `--snp` and `--matrix` take several files
and can be repeated. Logs go to stderr.

Option names follow ebscore: `--ponmatrix`, `-o/--outdir`, `--chrom`, `--bed`, `--pileup`, `-t`, `--tmpdir`, `--quiet`.
`--tumor-name/--normal-name` are labels for file names (ebscore's `--tumor/--normal` are BAM paths). Without them the
labels come from a legacy sample ID `<patient>_<tumor>-<normal>`; sarek pair IDs (`*_vs_*`) need explicit names.

## Failures instead of empty outputs

A task that cannot produce a meaningful result exits non-zero (2 for data problems, 1 for a failing shell pipeline)
with one log line that names the cause. Examples: a missing or mistyped `--gc`/`--genmap`, a track without data for a
chromosome of the sample, an empty pileup, a BED that is unsorted or does not match the pileup, a sample binned with
another BED or coordinate mode than the CNV PON (checked with `CNV/pon.json` and on the bins themselves), fewer than
90 % of the sample bins found in the PON (`--min-pon-match`), PON chromosomes missing from the sample (a partial
`combine`; `--allow-missing-chroms`), or no SNP left for ASCAT (`--allow-empty` writes header-only tables instead).

## Legacy compatibility

- `--legacy-coords` (`pon-coverage`, `tn2cnv`, `combine`, `sample`) reproduces the old bin labels and ExonPos. Use it
  only to combine new samples with PON coverage files built by the legacy tools; otherwise rebuild the CNV PON from
  `matrix/`. The dev PON `somvar_HAEv7_NovaSeq_v1` has legacy `cov/` and `CNV/` tables built on the **padded**
  SS_HAEv7 BED: run samples against it with `--bed SS_HAEv7_hg38_Padded.bed --legacy-coords`.
- `pon-combine --legacy-stats` reproduces the old PON tables (filter stats not recomputed).
- All deliberate changes and the measured deltas on a study sample: `CHANGELOG.md`.

## Layout

`src/mycnv/`: `raw` (tn2cnv/pon2cnv), `bed` (BED reading, bin ExonPos), `tracks` (GC/genmap reading), `trackgen` (GC/genmap building), `combine` (PON and
sample combine), `ponmeta` (PON provenance), `rolling`, `rolling_cov`, `rolling_snp`, `convert` (ASCAT input),
`pipeline` (steps used by the CLI), `cli`, `plot*`, `cluster` (GMM-caller groundwork; needs scikit-learn), `mawk/`.
Notebooks in `nb/` are the development history and still use the legacy `code/` imports (plan 2.3 moves them to
`import mycnv`). Tests: `tests/README.md`.
