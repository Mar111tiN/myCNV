# Changelog

## 2.0.0.dev0 (`kumodev`, 2026-11-22)

First packaged version, for the somvar → kumo-variant-calling-nf port (plan 3a). `main` (`db0801d`, the code the
somVarWES study runs used) stays the legacy reference. Parity tests compare with outputs of the unmodified `main` code
in the `somvar-legacy` env (`tests/legacy_outputs.py`, outputs in `tests/data/legacy/`).

### Structure
- Package `src/mycnv/` with the CLI `mycnv` (`pon-coverage`, `pon-combine`, `tn2cnv`, `combine`, `rolling`,
  `ascat-input`, `sample`, `config`). Parameters come from the packaged defaults (= somVarWES
  `configs/general/config_CNV.yaml`, `GCnormalize: True`), an optional `--config` YAML and CLI options.
- mawk tools are package data (`src/mycnv/mawk/`), run with `bash`: `cleanpileup`, `filterBed`, `pile2CNV`, `PON2CNV`.
- Removed on `kumodev`: `bamCoverage`, `rollingCoverage`, `cleanSNP`, `snpVAF`, `test` (2020 approach, unused), the
  `CNV_raw` LEGACY functions `get_heteroSNP`/`bam2coverage`, `snakemake/`, `env/`, `envMac/`, `.vscode/`.
- GC and mappability tracks are explicit inputs: a directory of per-chromosome files or one full-genome file.
- Per-sample steps run on the whole sample; `--chrom` restricts a run to one chromosome. `pon-coverage` runs one job
  per chromosome (`--threads` in parallel) (plan D16).

### Behaviour changes vs `main` (deliberate)
| Area | `main` | now | effect / how to get `main` |
|---|---|---|---|
| mawk shebang | `#!/bin/sh` with bash syntax: under dash all options were silently ignored | `#!/usr/bin/env bash`, always run with bash | — |
| `cleanpileup -s` | compared as string (ignored with ≥ 10 samples) | numeric | — |
| `filterBed` ExonPos | `BEDSUM − BEDSTART + pos`: not monotonic (NHL1 chr1: 8,482 decreasing steps), FullExonPos ties between regions | cumulative covered bases, contiguous | `-L` / `--legacy-coords` |
| `filterBed` | early `exit` after the last BED region; chromosome order had to match | any order; no early exit; BED header lines skipped | — |
| `pile2CNV`/`PON2CNV` bin labels | `Pos` = start of the *next* bin met; last bin of a chromosome its own start (duplicate label, 19 pairs in NHL1) | `Pos` = own bin start; bin ExonPos = ExonPos of its first covered base | `-L` / `--legacy-coords` |
| multi-chromosome streams | last bin of a chromosome printed under the next chromosome's name | flushed per chromosome | (main ran per chromosome) |
| `PON2CNV` first base of a bin | lost (undefined `$COL[s]`), ≈ 1/binSize of the coverage | counted | measured exactly in `tests/test_mawk.py` |
| PON outlier removal | stats not recomputed (filter == full stats); full tables got the outlier NaNs | filter stats recomputed on outlier-free coverage; full tables unmodified | `pon-combine --legacy` |
| FullExonPos offsets | lexicographic chromosome order (chr1, chr10, …, chr2) | natural order | — (rolling is per chromosome, so values do not depend on the order) |
| rolling windows | crossed chromosome boundaries (also the interpolation) | per chromosome; chromosomes shorter than the window get partial windows | — |
| ties in FullExonPos (legacy coords only) | interpolation neighbour chosen by numpy's unstable sort | (Chr, FullExonPos) stable sort, genomic neighbour | — |
| `diff_exp` | never read (code read `diffexp`, exponent always 2) | read (`diffexp` still accepted); study config: 4 → only `*Diff` columns change | set `diff_exp: 2` |
| `rolling().agg(ddof=)` | passed to mean/sum (TypeError on pandas ≥ 2) | only std/var | — |
| PON SNP merge | merged on all shared columns, incl. the float `map*` columns (breaks when files were parsed by different pandas versions) | on Chr, Pos, ExonPos | — |
| `Chr` categorical | `df.loc[:, "Chr"] = Categorical` (object dtype on pandas ≥ 2) | `df["Chr"] = …` | — |
| ASCAT names | `<patient>_<tumor>-<normal>` parsed on the first `_` | explicit `--tumor/--normal`, else the last `_` field | — |
| plots | `plt.cm.get_cmap`, `Series[-1]`, 23 rectangles vs fewer colours | `matplotlib.colormaps`, `.iloc`, one rectangle per present chromosome; empty tracks skipped | — |
| missing GC/genmap path | NameError | warning, empty columns | — |

### Measured on NHL1_A-B (whole sample, `raw/` → `data/*.cnv.snp.gz`, study config)
- `main` in `somvar-legacy` reproduces the study file exactly (12,939 rows, all 27 columns).
- The new code keeps the same 12,939 SNPs and identical values in all non-interpolated columns. SNPdensity/fallSNP
  shift by ≤ 0.0054 / 0.0028 (global min-max no longer sees cross-chromosome windows); SNPdensityDiff changes with
  `diff_exp` 4; the interpolated coverage columns (`log2ratio*`, `PONcov_*`, `GCratio`) differ on ~300–375 of 12,939
  rows (FullExonPos ties under the legacy coordinates and chromosome boundaries).
