# mycnv

PON-based copy-number profiles from a paired (normal, tumor) pileup: 10 bp coverage bins as log2ratio vs a Panel of
Normals, het-SNP B-allele frequencies, GC/mappability filters, rolling windows, and ASCAT input. This is the CNV core of
the somVarWES pipeline, packaged for the kumo-variant-calling-nf port (branch `kumodev`; `main` is the legacy code).

## Install

```bash
conda env create -f environment.yml -n mycnv        # unpinned conda-forge/bioconda env (mawk, pandas, matplotlib, ...)
conda run -n mycnv python -m pip install --no-deps --no-build-isolation .
conda run -n mycnv python -m pytest -q              # tests incl. parity with the legacy code
```
Container: `docker build -t mycnv:<ver> .` (named stages `build` → `test` → `runtime`; `--target test` runs the tests).

## Inputs

- **Paired pileup** (normal = sample 1): `samtools mpileup -f ref.fa -q 20 -Q 25 -l <covered BED> NORMAL TUMOR |
  cleanpileup -d` (as written by `ebscore pileup`), plain or gzipped. `--raw-pileup` accepts raw mpileup text.
- **PONmatrix** `matrix/{chrom}.pon.gz` (`Chr Start Ref A G C T I D Depth`, cells `s1|…|sN=s1|…|sN`), the pile2count
  matrix shared with ebscore.
- **BED**: the capture BED; use the same flavour for the PON and the samples (plan L9: `covered`).
- **GC** (`--gc`): directory with `{chrom}.gc100-10.gz` (`Start GC/AT`) or one file with a `Chr` column.
- **Mappability** (`--genmap`): directory with `hg38_genmap.{chrom}.txt.gz` (`Chr Pos 30_0 …`) or one file.

## Commands

```bash
# PON (once per kit + normal set): coverage/SNPs per chromosome, then the combined PON tables
mycnv pon-coverage --matrix-dir PONmatrix/matrix --bed covered.bed --gc gc100-10/ --genmap genmap/split/ \
    --out-dir PONmatrix --threads 8                     # -> PONmatrix/cov/{chrom}.cov.gz, snp/{chrom}.snp.gz
mycnv pon-combine --pon-dir PONmatrix                   # -> PONmatrix/CNV/pon.cov.{full,fullGC,filter,filterGC}.gz, pon.snp.gz

# per pair, all steps in one task
mycnv sample --pileup pair.pileup.gz --bed covered.bed --gc gc100-10/ --genmap genmap/split/ --pon-dir PONmatrix \
    --sample P1_T-N --out-dir CNV/P1 [--tumor T --normal N] [--plots]
# -> CNV/P1/raw/P1_T-N.{cov.gz,snp}  data/P1_T-N.{cnv.snp.gz,cnv.gz}  pre/P1_T-N_{N,T,germline}_{baf,logr}.tsv  [plots/]

# or step by step (same outputs)
mycnv tn2cnv   --pileup pair.pileup.gz --bed covered.bed --gc … --genmap … --out-cov s.cov.gz --out-snp s.snp [--chrom chr7]
mycnv combine  --cov s.cov.gz [more per-chromosome files] --snp s.snp [...] --pon-dir PONmatrix --out-cov raw.cov.gz --out-snp raw.snp
mycnv rolling  --cov raw.cov.gz --snp raw.snp --out-snp data.cnv.snp.gz --out-cnv data.cnv.gz
mycnv ascat-input --cnv data.cnv.snp.gz --sample P1_T-N --out-dir pre/
mycnv config   # print the effective parameters (defaults = the somVarWES study config)
```
Options shared by most commands: `--config my.yaml` (same structure as `mycnv config`), `-v` (debug logging, prints
the shell pipelines). Logs go to stderr.

## Legacy compatibility

- `--legacy-coords` (`pon-coverage`, `tn2cnv`, `sample`) reproduces the old bin labels and ExonPos. Use it only to combine
  new samples with PON coverage files built by the legacy tools; otherwise rebuild the CNV PON from `matrix/`.
- `pon-combine --legacy` reproduces the old PON tables (filter stats not recomputed).
- All deliberate changes and the measured deltas on a study sample: `CHANGELOG.md`.

## Layout

`src/mycnv/`: `raw` (tn2cnv/pon2cnv), `tracks` (GC/genmap), `combine` (PON and sample combine), `rolling`,
`rolling_cov`, `rolling_snp`, `convert` (ASCAT input), `pipeline` (steps used by the CLI), `cli`, `plot*`, `cluster`
(GMM-caller groundwork; needs scikit-learn), `mawk/`. Notebooks in `nb/` are the development history and still use
the legacy `code/` imports (plan 2.3 moves them to `import mycnv`).
