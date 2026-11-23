# Changelog

## 2.0.0.dev0 (`kumodev`, 2026-11-22)

First packaged version, for the somvar → kumo-variant-calling-nf port (plan 3a). `main` (`db0801d`, the code the
somVarWES study runs used) stays the legacy reference. Parity tests compare with outputs of the unmodified `main` code
in the `somvar-legacy` env (`tests/legacy_outputs.py`): public ones on synthetic inputs in `tests/data/legacy/`,
restricted ones on patient-derived slices in the git-ignored `tests/data/restricted/legacy/`.

### Structure
- Package `src/mycnv/` with the CLI `mycnv` (`pon-coverage`, `pon-combine`, `tn2cnv`, `combine`, `rolling`,
  `ascat-input`, `sample`, `config`). Parameters come from the packaged defaults (= somVarWES
  `configs/general/config_CNV.yaml`, `GCnormalize: True`), an optional `--config` YAML and CLI options.
- mawk tools are package data (`src/mycnv/mawk/`), run with `bash`: `cleanpileup`, `filterBed`, `pile2CNV`, `PON2CNV`.
- Removed on `kumodev`: `bamCoverage`, `rollingCoverage`, `cleanSNP`, `snpVAF`, `test` (2020 approach, unused), the
  `CNV_raw` LEGACY functions `get_heteroSNP`/`bam2coverage`, `snakemake/`, `env/`, `envMac/`, `.vscode/`.
- GC and mappability tracks are explicit inputs: a directory of per-chromosome files or one file.
- New `mycnv tracks gc|genmap` (`mycnv.trackgen`, plan D22): rebuilds both tracks from a FASTA for genome-bundle
  builds, because the studies' kit-specific tracks are lost. GC reproduces the old static track wherever the FASTAs
  agree (chr21: 99.985 % of rows, same format). Mappability uses genmap with the old column layout and "blank =
  unchanged" encoding; the existing reader decodes it unchanged. genmap and samtools are now in `environment.yml`.
  Provenance of the old generic genmap track (chr21-only index; a genome-wide check is not allowed on the
  workbench): 66–76 % of positions equal per mode, and 97.6 % agree where the old track says unique. This is
  consistent with genmap defaults on an exome-wide selection, but on a genome-wide index of a FASTA other than the
  GATK bundle (the bundle masks ~2 Mb of chr21 that the old GC track has values for).
- Per-sample steps run on the whole sample. `tn2cnv --chrom` (repeatable) is the optional split; its outputs all go
  into one `combine`. `sample` has no `--chrom`: GC normalisation and the log2ratio levels need the whole sample.
  `pon-coverage` runs one job per chromosome (`--threads` in parallel) (plan D16).
- CLI vocabulary shared with ebscore: `--ponmatrix DIR`, `-o/--outdir` (required wherever a folder is written: inputs
  are never modified), repeatable `--chrom`, `-t/--threads`, `--tmpdir`, `-v`, `--quiet`; `--tumor-name/--normal-name`
  for the labels. No argparse prefix matching. Data errors exit 2, failing shell pipelines exit 1.
- The image no longer needs seaborn (and scipy): the plots use matplotlib's `seaborn-v0_8` style; matplotlib is only
  imported for plots. `default_config.yaml` lost the unused `pileup` and `ascat` sections; `--config` reports and ignores
  keys mycnv does not read.

### Behaviour changes vs `main` (deliberate)
| Area | `main` | now | effect / how to get `main` |
|---|---|---|---|
| mawk shebang | `#!/bin/sh` with bash syntax: under dash all options were silently ignored | `#!/usr/bin/env bash`, always run with bash | — |
| `cleanpileup -s` | compared as string (ignored with ≥ 10 samples) | numeric | — |
| `filterBed` ExonPos | `BEDSUM − BEDSTART + pos`: not monotonic (NHL1 chr1: 8,482 decreasing steps), FullExonPos ties between regions | cumulative covered bases, contiguous | `-L` / `--legacy-coords` |
| `filterBed` | early `exit` after the last BED region; chromosome order had to match | any order; no early exit; BED header lines skipped | — |
| `filterBed` BED order | an unsorted chromosome or one that came back after another one dropped or misassigned positions silently | error (`sort -k1,1 -k2,2n`) | — |
| overlapping BED regions | every region added `end − start + 1` to ExonPos, so overlapping (and bookended, which share a base under the start..end-inclusive convention) bases were counted twice: vendor padded BEDs are not merged (SS_HAEv7 padded: 29,488 overlaps, 115 bookended; the tier-A `testkit.padded.bed`: 12) | fixed coordinates merge overlapping/touching regions per chromosome before summing (filterBed v2.2.0, `bed.merge_regions`, also used for the bin ExonPos and the BED digests); a BED and its merged copy give identical outputs | `-L` / `--legacy-coords` keeps the unmerged regions |
| `pile2CNV`/`PON2CNV` bin labels | `Pos` = start of the *next* bin met; last bin of a chromosome its own start (duplicate label, 19 pairs in NHL1) | `Pos` = own bin start | `-L` / `--legacy-coords` |
| bin ExonPos | first streamed base's ExonPos − pos % 10 | ExonPos of the bin's first **BED** base (`bed.bin_exon_pos`), the same for a sample and a PON whatever bases their streams hold | `--legacy-coords` |
| sample ↔ PON join | inner join on Chr, Pos, ExonPos, unchecked (0 bins and exit 0 when BED or coordinate mode differed) | on (Chr, Pos); every shared bin must have the same ExonPos; ≥ 90 % of the sample bins must find a PON bin; `CNV/pon.json` mode/BED check | legacy mode keeps the 3-column join |
| multi-chromosome streams | last bin of a chromosome printed under the next chromosome's name | flushed per chromosome | (main ran per chromosome) |
| `PON2CNV` first base of a bin | lost (undefined `$COL[s]`), ≈ 1/binSize of the coverage | counted | measured exactly in the tests |
| `pile2CNV` het-SNP alt reads | A/C/G/T, D/d and only forward-strand insertions `I`: het insertions got ~half their VAF | also reverse-strand `i` | — |
| PON outlier removal | stats not recomputed (filter == full stats); full tables got the outlier NaNs | filter stats recomputed on outlier-free coverage; full tables unmodified | `pon-combine --legacy-stats` |
| sample GC normalisation | over chr1–22, X (`combine_chrom_cnv` loaded only those) | the same: bins on chromosomes without PON coverage (chrY, alts) are dropped before normalising | — |
| FullExonPos offsets | lexicographic chromosome order (chr1, chr10, …, chr2) | natural order | — (rolling is per chromosome, so values do not depend on the order) |
| rolling windows | crossed chromosome boundaries (also the interpolation) | per chromosome; chromosomes with fewer than 2 · (window − 10) − 1 rows get partial windows (below that the full-window margins left a block of NaN) | — |
| ties in FullExonPos (legacy coords only) | interpolation neighbour chosen by numpy's unstable sort | (Chr, FullExonPos) stable sort, genomic neighbour | — |
| `diff_exp` | never read (code read `diffexp`, exponent always 2) | read (`diffexp` still accepted); study config: 4 → only `*Diff` columns change | set `diff_exp: 2` |
| `rolling().agg(ddof=)` | passed to mean/sum (TypeError on pandas ≥ 2) | only std/var | — |
| fallSNP without any PON VAF | 0/0: NaN, every SNP dropped | 0 (no SNP is a PON artefact) and a warning | — |
| PON SNP merge | merged on all shared columns, incl. the float `map*` columns (breaks when files were parsed by different pandas versions) | on Chr, Pos, ExonPos | — |
| `Chr` categorical | `df.loc[:, "Chr"] = Categorical` (object dtype on pandas ≥ 2) | `df["Chr"] = …` | — |
| header-only per-chromosome tables | (pandas 1.5 kept float dtypes) | skipped before concatenation (pandas 3 turned every column into object) | — |
| ASCAT names | `<patient>_<tumor>-<normal>` parsed on the first `_` | `--tumor-name/--normal-name` (both or neither), else the last `_` field of a legacy ID; `*_vs_*` IDs rejected | — |
| ASCAT input without SNPs or logR | header-only / all-NaN tables, exit 0 | error; `--allow-empty` writes them | — |
| plots | `plt.cm.get_cmap`, `Series[-1]`, 23 rectangles vs fewer colours; crash on a one-SNP track | `matplotlib.colormaps`, `.iloc`, one rectangle per present chromosome; empty tracks skipped; tick guard; written after the ASCAT input | — |
| missing / mistyped GC or genmap track | NameError | error (`--gc` required; `--genmap` required for samples); contigs outside chr1–22/X only warn; NaN-map SNPs are counted in the log | — |
| empty pileup, no bin in the BED | header-only outputs | error | — |

### Known differences that are not code changes
- **genmap track**: the studies used a kit-specific split (`HAEv7split`) that is not available. The generic
  `genmap/split/` track differs: on NHL1 chr7, 572 of 1,487 SNPs get other map values and 20 fall into gaps of the
  change-point-encoded track (no value, filtered out); with the default thresholds 1,233 SNPs pass with the study
  values and 1,289 with the generic track (100 study-only, 156 generic-only). Map-filter parity with the studies needs
  the kit-specific track or a dense per-base track. The gaps are decoded as before (exact position match after
  ffill/bfill of the stored rows); a `merge_asof` decoding matched the study values for only 5 of the 20 gap SNPs.
- **Diff scaling on multi-chromosome samples**: `Diff = |R − L| / max` still uses the sample-wide max, as in main.
  main's cross-boundary windows could set that max; with per-chromosome windows a sample whose chromosome levels jump
  gets a smaller max, so all smoothed values shift slightly (synthetic 3-chromosome fixture: most coverage rows;
  NHL1: see below). The SNP-density min-max normalisation likewise no longer sees main's near-zero cross-boundary
  density, so a few more SNPs can pass the density filter (NHL1 chr7+chr21 slice: 41 instead of 36, all 36 kept).

### Measured on NHL1_A-B (whole sample, `raw/` → `data/*.cnv.snp.gz`, study config)
- `main` in `somvar-legacy` reproduces the study file exactly (12,939 rows, all 27 columns).
- The new code keeps the same 12,939 SNPs and identical values in all non-interpolated columns. SNPdensity/fallSNP
  shift by ≤ 0.0054 / 0.0028 (global min-max no longer sees cross-chromosome windows); SNPdensityDiff changes with
  `diff_exp` 4; the interpolated coverage columns (`log2ratio*`, `PONcov_*`, `GCratio`) differ on ~300–375 of 12,939
  rows (FullExonPos ties under the legacy coordinates and chromosome boundaries).
