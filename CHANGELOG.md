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
- Coverage sex calls (`mycnv.sex`, design items M4/M5, plan N-Q7/D-Q3): ratio = `2 ** (median − median)` of chrX
  non-PAR vs autosomal bins of a log2ratio column against an effectively female reference (a plain coverage column:
  median(chrX non-PAR) / median(autosomes), kit-dependent, only the PON's start); NaN/inf and zero-coverage bins
  ignored; GRCh38 PARs (chrX:10,001–2,781,479, 155,701,383–156,030,895) excluded; XY < 0.65 < ambiguous < 0.85 < XX;
  unknown below 100 non-PAR chrX bins or without autosomal bins. New config section `sex:` (`xy_max`, `xx_min`,
  `min_bins`), checked when the config is loaded: a bad section, also one that is not a mapping (`sex: 0.9`), fails
  every command with exit 2 before it does work.
  - `pon-combine`: every PON normal is called on its mean-normalised coverage (before doubling) as a log2ratio vs the
    per-bin mean of the other normals with the XY normals' chrX doubled. The XY set is recomputed until it is stable.
    It starts from the split of the sorted log2 plain coverage ratios whose lower part, doubled, leaves the PON
    tightest (least sum of squares), which does not depend on the kit's chrX mix. A PON of one sex (no split beats
    doubling none) starts from the plain ratio cut at `(xy_max + xx_min) / 2`, with a warning. `CNV/pon.json` gains
    `sex` (per PON column: `column`, `name` from `PON_list.txt` if it matches the columns, `chrx_ratio`, `call`,
    `doubled`) and `sex_thresholds` (thresholds, PARs, `pon_call`, `doubling`).
  - `sample` and `combine` write `sex.tsv` (`sample role column chrx_ratio log2_diff n_chrx_bins n_auto_bins call
    ascat_gender`) from the raw coverage table; `ascat_gender` is the normal's call (ambiguous/unknown → XX, with a
    warning) and is the gender the pipeline's ASCAT step uses (D-Q3). `combine` writes it into the `--out-cov` folder
    (`--out-sex FILE` elsewhere); the new subcommand `mycnv sex --cov raw/<s>.cov.gz -o sex.tsv` makes the same file.
    An empty raw table, or one without a finite autosomal log2ratio, exits 2. The call is made before any output is
    written. The existing `raw/`, `data/`, `pre/` outputs are unchanged (byte-identical on the synthetic
    multi-chromosome kit, the synthetic CNV PON + sample and tier A).
- Tumor-only samples (design item M1, plan D27, N-Q3, D-Q3 mode B): `sample --tumor-only --snp-sites SITES` (and
  `tn2cnv --tumor-only --snp-sites`, `combine`/`rolling`/`sex`/`ascat-input --tumor-only` for the step-by-step path).
  main had no tumor-only CNV path.
  - Input: the one-sample pileup of `ebscore pileup --tumor` (`Chr Start Ref Depth1 Read1`, tumor = sample 1) and a
    sites file: VCF/VCF.gz (CHROM POS REF ALT; biallelic SNVs kept; indels and symbolic alleles dropped; a multi-allelic
    site dropped both as one record and split into biallelic records at one position, an identical record counted
    once) or a 2-column `Chr Pos` TSV (1-based; `#` lines are comments, a first line whose Pos has no digit is the
    header, a bad Pos is an error; a BED is rejected). VCF = a `##fileformat=VCF` line or a `#CHROM POS ID REF ALT`
    header. Read in chunks of 500,000 rows and only kept at candidate positions (`mycnv.sites`). The file is checked
    before the pileup pass (exists, regular file, at least one record).
  - Coverage bins: the same `filterBed | pile2CNV` pass and `bed.bin_exon_pos` as a pair (`raw._pile2cnv`, shared by
    `tn2cnv` and the new `tumor2cnv`): on the tumor columns of a pair the bins are identical, `Cov1` = the pair's
    `Cov2`.
  - Het-SNP candidates: the sites where the tumor has depth ≥ `hetSNP.minDepth` and a VAF within the new config key
    `hetSNP.tumorOnlyVAF` (default `[0.1, 0.9]`, inclusive; `--tumor-vaf`; owner-tunable, no legacy value), with
    pile2CNV's alt reads (A/C/G/T, D/d, I/i). `maxPONVAF`, mappability and density filters apply downstream as for
    pairs.
  - Downstream with one sample column (`Cov1`, `VAF1`, `log2ratio1` = tumor): combine, rolling and the ASCAT tables
    work per column, so a tumor-only run on a pair's tumor columns gives the pair's tumor values exactly (synthetic CNV
    PON + sample: raw, rolled and ASCAT tables identical). `write_ascat(tumor_only=True)` writes
    `pre/<S>_<T>_{baf,logr}.tsv` + the mock germline `pre/<S>_germline_{baf,logr}.tsv` (BAF 0.5, logR 0), no normal
    tables; the tumor label defaults to `--sample`. `sex.tsv` has one row (role `tumor`, `log2ratio1`) and
    `ascat_gender` = the tumor's call; an XY tumor-only call logs a warning (a female tumor with a chrX loss looks XY).
    Plots are written for the tumor only.
  - CLI rules (exit 2): `--tumor-only` needs `--snp-sites`; `--snp-sites`/`--tumor-vaf` need `--tumor-only`;
    `--normal-name`, `--normal-vaf` and `--no-mock-germline` are rejected with `--tumor-only` (D-Q3: one mock-germline
    path). A two-sample pileup with `--tumor-only` and tables of a pair in the tumor-only commands fail, and vice
    versa (see the table below); a sites file sharing no chromosome with the BED fails (naming).
  - The paired path is unchanged: its outputs are byte-identical on the synthetic multi-chromosome kit and the
    synthetic CNV PON + sample (parity tests) and on tier A (`raw/ data/ pre/` content-identical to `expected/`).
- Test fixtures: the synthetic CNV PON (`tests/data/synth/pon`, and the sample on its bins) had all chrX bins at
  100–110 kb, inside PAR1, where an XY normal has two copies; the PAR-aware call found 0 non-PAR bins there. Its chrX
  now starts at 3.0 Mb (`tests/synth.py` `SYN_START`; same random draws). The `main` outputs made from it
  (`tests/data/legacy/{pon_synth_CNV,combine_synth_*,rolling_synth_all*}`) were regenerated with unmodified `main` in
  `somvar-legacy` and are identical to the previous ones up to chrX `Pos` + 2,900,000.
- CLI contract snapshot (plan C1, 2026-11-24). This is a test only: the CLI and its outputs do not change.
  `tests/data/cli_contract.json` (schema `mycnv-cli-contract/1`, the layout of ebscore's C1 snapshot) is written by
  `tests/make_cli_contract.py` from `mycnv.cli.build_parser()`. It covers every subcommand, the nested `tracks gc` and
  `tracks genmap` included, with the sorted option strings. Per option it records the flags, required, the action kind,
  `nargs`, `type`, `choices` and the default when that is a constant. It also records the mutually exclusive groups,
  the positionals, and the `--version` format with the version masked (`mycnv <version>`). The sex (`sex`,
  `combine --out-sex`) and tumor-only (`--tumor-only`, `--snp-sites`, `--tumor-vaf`) options are in it. Help texts are
  not: argparse formats `--help` differently per Python version, and the snapshot is the same on Python 3.9 and 3.15.
  - `tests/test_cli_contract.py` fails on any difference, compared as canonical JSON (so `1`, `1.0` and `true`
    differ). The message lists each changed entry and the kumo-variant-calling-nf modules that call the changed
    subcommand (`MYCNV_PONCOVERAGE`, `MYCNV_PONCOMBINE`, `MYCNV_SAMPLE`; design §3). It says to rerun
    `python tests/make_cli_contract.py` and review those modules. When the imported `mycnv` is a stale install (another
    copy than `src/`, with other code), it says that instead.
  - The test also parses the design §3 command lines (`pon-coverage`, `pon-combine`, `sample` paired and
    `--tumor-only`) with prefix matching off and runs the CLI's cross-option checks on them. It checks that no
    (sub)command parser allows prefix matching (plan risk "moving CLIs").

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
| PON chrX doubling (`amazonize`) | chrX ×2 (all chrX bins, PAR included) for normals whose normalised chrX mean (PAR included) is < 75, decided separately in the GC and the mean table; not recorded. Fails with column dropouts (HAEv7, chr21 + chr22 + chrX: 8 normals in the GC table, none in the mean table; the coverage call has the same 15 XY on every subset) | the coverage sex call of each normal decides for both tables: log2ratio vs the per-bin mean of the other normals (XY chrX doubled, iterated; PAR and zero bins excluded), the thresholds of the samples. XY: only the chrX non-PAR bins doubled (two PAR copies in XY; main put the PON PAR ~0.65 log2 above the rest of chrX, an apparent loss over PAR1 in every sample); XX, ambiguous, unknown not doubled (warning per ambiguous/unknown normal, a summary warning, and a warning when main's rule differs); recorded in `CNV/pon.json`. Same columns as main's rule on the synthetic PON (Cov3, Cov6), on the restricted `pon_mini` (none: its chrX slice lies in PAR1) and on the HAEv7 PON (15 of 28, see below) | `pon-combine --legacy-stats` (main's rule, all chrX bins; the calls are only recorded) |
| sample GC normalisation | over chr1–22, X (`combine_chrom_cnv` loaded only those) | the same: bins on chromosomes without PON coverage (chrY, alts) are dropped before normalising | — |
| GC class with autosomal mean 0 | factor 100/0: its autosomal bins NaN (0 · inf), its chrX bins with coverage **inf**, which spread into the PON stats and rolling windows (seen when the PON was combined on chr7, chr21, chrX only: 15 inf values) | factor NaN: all its bins NaN, with a warning; the bins are kept | — (whole-genome PON and samples: no such class in HAEv7) |
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
| pileup layout | not checked: a one-sample pileup ran as a "pair" (all downstream steps on one column) and a raw mpileup read as a clean one gave garbage | the first line is checked: a clean pileup needs the `cleanpileup -d` header, a raw one (`--raw-pileup`) 3 + 3 per sample columns; a pair needs 2 samples, `--tumor-only` 1; the error names the flag. A regular file is checked before the pass; a pipe or FIFO (`/dev/stdin`), which can be read only once, on the pile2CNV stream header after it (same errors, exit 2; a clean pileup on a pipe with `--raw-pileup` is not detected). A missing pileup is now exit 2 (data problem) instead of 1 (failing `cat`) | — (valid pairs, also streamed: unchanged outputs) |
| sample tables of `combine`/`rolling` | any number of sample columns | a pair needs `Cov1`/`Cov2` (and `VAF1`/`VAF2`, `log2ratio1`/`log2ratio2` where present); `--tumor-only` exactly `*1` | — |

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

### Measured on AMLMono7 03_A-B (notebooks 3 → 4 → 5 on the study `raw/`, 2026-11-23)
- `main` in `somvar-legacy` reproduces the study `data/03_A-B.cnv.snp.gz` exactly (11,945 rows, all 27 columns).
- The new code keeps 11,995 SNPs: 51 more and 1 fewer. Every `fallSNP` is ≈ 0.575 × main's (p01–p99 0.533–0.580),
  because `fallSNP = SNPdensity · PONVAF / max` is scaled by one SNP: chr7:100,958,115 in MUC3A, where two SNPs 625 bp
  apart share FullExonPos 43,183,799 under the legacy coordinates. The stable (Chr, FullExonPos) sort orders such ties
  differently, the 5-SNP density window there changes, and the genome-wide max rises 1.74×. SNPs with main's fallSNP
  just above `maxFallSNP` 0.01 then pass. Fixed coordinates have no such ties; the max-scaling itself (one SNP cluster
  sets the cut for the whole sample) is unchanged legacy behaviour.
- Rolled coverage: median |Δ| ≈ 1e-7, 289 rows > 0.01 (as on NHL1).

### Sex calls measured on the HAEv7 PON `somvar_HAEv7_NovaSeq_v1` (D21 `cov/`, 28 normals; 2026-11-23, not genome-wide)
- main's rule doubles 15 of the 28 normals, the same 15 in the GC and the mean table, whether chr7 + chr21 + chrX or
  chr1 + chr2 + chr7 + chr21 + chrX are combined.
- The plain coverage ratio of a PON column (no per-bin reference) is not centred on 1 and 0.5 for this kit: XX
  1.15–1.35, XY 0.60–0.68 (quotient ≈ 0.5); zero-coverage bins raised it further (chr22 has column dropouts: 10 of
  the 28 normals have 29–96 % zero bins there). Called on that ratio with the default thresholds, chr7 + chr21 + chrX
  gave 6 XY, 9 ambiguous (not doubled), 13 XX, and chr21 + chr22 + chrX left 10 normals unknown (first version of this
  change, review 2026-11-23).
- Called as now (log2ratio vs the other normals' per-bin mean, XY chrX doubled; zeros excluded), with the default
  thresholds:
  - chr7 + chr21 + chrX: 15 XY (0.483–0.526), 13 XX (0.921–1.032), none ambiguous; the 15 doubled normals are main's
    15. After doubling, the PON PAR sits at the level of the rest of chrX (fullGC, PON
    mean vs the mean of the XX normals, median log2: PAR +0.013, chrX non-PAR +0.020, autosomes 0.000); main's
    doubling of the PAR put it about 0.65 log2 above the rest of chrX.
  - chr21 + chr22 + chrX: the same 15 XY (0.396–0.524), 9 XX, 4 XX-level normals ambiguous (0.73–0.85, not doubled,
    with warnings); main's rule doubles 8 normals in the GC table and none in the mean table here.
  - The start matters only at the edges: in an exploratory run (all-column mean as reference), the starts "all",
    main's rule and the plain-ratio cut reached the same 15 XY on both subsets, while an empty start stayed empty on
    chr7 + chr21 + chrX (against a mixed reference the XY normals fall into the ambiguous band). The reference now
    leaves the called normal out, so its call does not depend on its own doubling (with N = 8 that dependence made a
    second fixed point in the tests). Genome-wide values were not computed on the workbench (memory cap).
  - Start rule (review 2026-11-24). The first version started from the plain-ratio cut alone. That fails on a balanced
    mixed PON (N-Q7's IMMUcan case, 50 XX + 50 XY) on a kit with a high chrX mix, where the plain XY ratio is above
    about 0.76. Nothing is seeded, the XY normals sit at about 0.5 / 0.75 ≈ 0.67 against the undoubled mixed reference
    (ambiguous), and that empty set is a fixed point. Every XY normal stayed undoubled, with warnings only. Synthetic
    PON, 28 normals / 15 XY and 100 / 50, chrX mix 1.6 or 1.8: 0 XY, all XY ambiguous.
    - The fix is the split start. A second start from "every normal XY", keeping the fixed point with fewer ambiguous
      normals, was rejected. With one genuinely ambiguous normal in an XX-majority PON, "every normal XY" is a fixed
      point with no ambiguous normal, so every normal would be doubled.
    - Synthetic PONs (28 / 15, 100 / 50, 28 / 6, 8 / 2, 8 / 4, 28 / 1, 28 / 27, chrX mix 0.8–1.8, also with one normal
      at 75 % chrX): the split start gives the true XY set every time.
    - The remaining limit is a PON of one sex: all XY on a kit with a plain XY ratio above 0.75 is called XX.
    - HAEv7: the split start already is the final XY set, and the calls on both subsets are the same as from the cut.
      The restricted `pon_mini` is all unknown (no non-PAR chrX bins), as before.
