# mycnv

PON-based copy-number profiles from a paired (normal, tumor) pileup, or from a tumor alone (`--tumor-only`, het SNPs
at given common-SNP sites): 10 bp coverage bins as log2ratio vs a Panel of
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
with a contract check, the tier-A fixture smoke inside the image (`tests/smoke/smoke_mycnv.sh` on
`tests/data/plumbing/`, as the image user and as an arbitrary uid), and Trivy 0.69.3: a Dockerfile config scan, the vulnerability gate (library packages, fixed
vulnerabilities, HIGH/CRITICAL; accepted findings go into `.trivyignore` with a reason) and an informational report
that includes the OS packages.

## Inputs

- **Paired pileup** (normal = sample 1): `samtools mpileup -f ref.fa -q 20 -Q 25 -l <padded BED> NORMAL TUMOR |
  cleanpileup -d` (as written by `ebscore pileup`), plain or gzipped. `--raw-pileup` accepts raw mpileup text. The
  first line is checked: a pair needs 2 samples; a one-sample (tumor) pileup needs `--tumor-only` (next section).
  A regular file is checked before the pass; a pipe or FIFO (`--pileup /dev/stdin`) can be read only once, so it is
  checked on the stream header after the pass (same errors). One gap remains there: a clean pileup on a pipe given
  with `--raw-pileup` is not detected, so give `--raw-pileup` only for raw text.
- **Tumor-only** (`--tumor-only`): the one-sample pileup of `ebscore pileup --tumor T` (`Chr Start Ref Depth1 Read1`,
  tumor = sample 1) and a common-SNP sites file (`--snp-sites`), see "Tumor-only samples".
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
# -> CNV/P1/raw/P1_T-N.{cov.gz,snp}  sex.tsv  data/P1_T-N.{cnv.snp.gz,cnv.gz}  pre/P1_T-N_{N,T,germline}_{baf,logr}.tsv
#    [plots/]

# or step by step (same outputs; combine writes sex.tsv into the --out-cov folder unless --out-sex says otherwise);
# a per-chromosome split stops after tn2cnv, `combine` takes all its outputs
mycnv tn2cnv   --pileup pair.pileup.gz --bed kit.padded.bed --gc … --genmap … --out-cov s.chr7.cov.gz --out-snp s.chr7.snp \
    [--chrom chr7 ...]
mycnv combine  --cov s.chr*.cov.gz --snp s.chr*.snp --ponmatrix ponmatrix_cnv --out-cov raw.cov.gz --out-snp raw.snp \
    [--out-sex sex.tsv] [--tumor-name T --normal-name N]                                   # -> + sex.tsv
mycnv sex      --cov raw.cov.gz -o sex.tsv [--tumor-name T --normal-name N | --sample P1_T-N]   # = combine's sex.tsv
mycnv rolling  --cov raw.cov.gz --snp raw.snp --out-snp data.cnv.snp.gz --out-cnv data.cnv.gz
mycnv ascat-input --cnv data.cnv.snp.gz --sample P1_T-N -o pre/
mycnv config   # print the effective parameters (defaults = the somVarWES study config)

# tumor-only sample (no normal): one-sample pileup + common-SNP sites; the tumor label defaults to --sample
mycnv sample --tumor-only --snp-sites kit.snp_sites.vcf.gz --pileup T1.pileup.gz --bed kit.padded.bed --gc gc100-10/ \
    --genmap genmap/split/ --ponmatrix ponmatrix_cnv --sample T1 -o CNV/T1 [--tumor-name T] [--plots]
# -> CNV/T1/raw/T1.{cov.gz,snp}  sex.tsv (one row)  data/T1.{cnv.snp.gz,cnv.gz}
#    pre/T1_T1_{baf,logr}.tsv + pre/T1_germline_{baf,logr}.tsv (mock germline)  [plots/]
# step by step: every per-sample command takes --tumor-only (tn2cnv also --snp-sites)
mycnv tn2cnv --tumor-only --snp-sites kit.snp_sites.vcf.gz --pileup T1.pileup.gz --bed … --gc … --genmap … \
    --out-cov t.cov.gz --out-snp t.snp [--chrom chr7 ...]
mycnv combine --tumor-only --cov t*.cov.gz --snp t*.snp --ponmatrix ponmatrix_cnv --out-cov raw.cov.gz --out-snp raw.snp \
    --sample T1
mycnv rolling --tumor-only --cov raw.cov.gz --snp raw.snp --out-snp data.cnv.snp.gz --out-cnv data.cnv.gz
mycnv ascat-input --tumor-only --cnv data.cnv.snp.gz --sample T1 -o pre/
mycnv sex --tumor-only --cov raw.cov.gz -o sex.tsv --sample T1
```
Every command takes `--config my.yaml` (structure of `mycnv config`; unknown keys are reported and ignored),
`-t/--threads` (parallel only in `pon-coverage`), `-v/--verbose` (prints the shell pipelines) and `--quiet`. Commands
that run the mawk tools take `--tmpdir`. `--chrom` is repeatable; `--cov`, `--snp` and `--matrix` take several files
and can be repeated. Logs go to stderr.

Option names follow ebscore: `--ponmatrix`, `-o/--outdir`, `--chrom`, `--bed`, `--pileup`, `-t`, `--tmpdir`, `--quiet`.
`--tumor-name/--normal-name` are labels for file names (ebscore's `--tumor/--normal` are BAM paths). Without them the
labels come from a legacy sample ID `<patient>_<tumor>-<normal>`; sarek pair IDs (`*_vs_*`) need explicit names.

## Tumor-only samples (`--tumor-only`, plan D27)

A sample without a normal runs through the same steps with one sample column, the tumor (`Cov1`, `VAF1`,
`log2ratio1`), against the same CNV PON. Only the het-SNP selection differs, because no normal says which positions
are germline het SNPs (N-Q3):
- **Coverage bins** come from the same `filterBed | pile2CNV` pass as a pair's (the tumor's `Cov1` equals a pair's
  `Cov2` on the same reads), with the same BED-derived ExonPos. The tumor logR is the log2ratio vs the PON mean
  (D-Q3).
- **Het-SNP candidates** are the positions of `--snp-sites` where the tumor has depth ≥ `hetSNP.minDepth` (15) and a
  VAF within `hetSNP.tumorOnlyVAF` (default `[0.1, 0.9]`, inclusive; `--tumor-vaf MIN MAX` or `--config`). Alt reads
  are counted as for pairs (A/C/G/T, D/d, I/i). The band is a new, owner-tunable parameter: a wider band keeps more
  LOH-shifted sites and more artefacts. The PON artefact filter (`filter.snp.maxPONVAF` on `CNV/pon.snp.gz`), the
  mappability and the SNP-density filters apply downstream as for pairs.
- **`--snp-sites`**: a VCF or VCF.gz, or a TSV with the two columns `Chr Pos` (1-based; a BED is rejected), plain
  or gzipped/bgzipped. The format is detected from the content: a file with a `##fileformat=VCF` line or a
  `#CHROM POS ID REF ALT …` header is a VCF.
  - VCF: only CHROM, POS, REF and ALT are read. Biallelic SNVs are kept; indels and symbolic alleles are dropped.
    A multi-allelic site is dropped whether it is one record (`ALT` `A,T`) or split into biallelic records at one
    position (`bcftools norm -m-`, gnomAD v4): a position with records of more than one allele is dropped, an
    identical record listed twice counts once.
  - TSV: `#` lines are comments, so the `#CHROM POS` of `grep -v '^##' sites.vcf | cut -f1,2` is fine. A first line
    whose Pos has no digit (`Chr Pos`) is the header. Every other line is data, and a Pos that is not a 1-based
    integer is an error.
  - The file is checked before the pileup pass: it must exist, be a regular file (it is read more than once: no
    pipe or process substitution) and hold at least one record. Chromosome names must be written as in the pileup and
    the BED (no `chr` prefix is added); a sites file that shares no chromosome with the BED fails.
  - The pipeline's file is the kit's `<kit>.snp_sites.vcf.gz`: common sites of gnomAD and dbSNP (or the one
    available), cut to the kit's padded BED (design K1). The file is read in chunks, so a large one is fine.
- **Outputs**: `pre/<S>_<T>_{baf,logr}.tsv` and the mock germline `pre/<S>_germline_{baf,logr}.tsv` (BAF 0.5, logR 0:
  ASCAT mode B, D-Q3), no normal tables; `<T>` is `--tumor-name`, default `--sample`. `sex.tsv` has one row (role
  `tumor`, column `log2ratio1`), and its `ascat_gender` is the tumor's own call. That call has a lower confidence than a
  pair's: a female tumor that lost one X has the chrX level of an XY sample and is called XY (logged as a warning; the
  pipeline cross-checks the sheet sex). Plots (`--plots`) are written for the tumor only.
- **Checks** (exit 2): `--tumor-only` needs `--snp-sites`, and `--snp-sites`/`--tumor-vaf` need `--tumor-only`;
  `--normal-name`, `--normal-vaf` and `--no-mock-germline` are rejected with `--tumor-only`. A two-sample pileup with
  `--tumor-only`, a one-sample pileup without it, and tables of the other layout in `combine`/`rolling`/`sex`/
  `ascat-input` fail with a message that names the flag.

## Sex calls (chrX of the PON normals and the samples)

One coverage-based call (`mycnv.sex`) is used for the PON and for the samples (plan N-Q7, D-Q3):
- `ratio = 2 ** (median(chrX non-PAR) − median(autosomes))` of a log2ratio column against an effectively female
  reference: XX ≈ 1, XY ≈ 0.5. (A plain coverage column gives `median(chrX non-PAR bins) / median(autosomal bins)` on the
  same scale, but without a per-bin reference it carries the kit's chrX/autosome target mix; the PON uses it only for
  the start, see below.) NaN/inf values and zero-coverage bins are ignored (a log2ratio is NaN where the coverage is 0,
  and a coverage column drops its 0 bins the same way). The PARs are GRCh38's (chrX:10,001–2,781,479 and
  155,701,383–156,030,895, 1-based inclusive), matched on the bin start `Pos`; autosomes are chr1–chr22.
- Call: **XY** if ratio < `sex.xy_max` (0.65), **XX** if ratio > `sex.xx_min` (0.85), else **ambiguous**; **unknown**
  with fewer than `sex.min_bins` (100) non-PAR chrX bins or no autosomal bins (e.g. a chr21-only test kit). Change the
  thresholds with `--config` (section `sex:`, a mapping, see `mycnv config`); a bad section fails every command
  before it starts (exit 2).

**PON** (`pon-combine`): every normal is called like a sample, on the log2ratio of its mean-normalised coverage (before
any doubling) against the per-bin mean of the *other* normals, in which the XY normals have their chrX doubled. The XY
set is recomputed until it no longer changes (usually one round). It starts from the plain coverage ratios. The kit's
chrX target mix scales all of them alike, and the XY normals of a mixed PON sit a factor 2 below the rest. So the start
is the split of the sorted log2 ratios whose lower part, doubled, leaves the PON tightest. This works for any kit and
any XX/XY mix. If no split beats doubling none, the PON is all one sex: the ratios alone cannot tell XX from XY, and the
start is the plain ratio cut at `(xy_max + xx_min) / 2` (with a warning). That cut assumes a chrX/autosome mix
near 1. On a kit whose plain XY ratio is above 0.75, an all-XY PON is called XX and not doubled; a single XX normal
among the XY ones avoids that. The chrX **non-PAR** bins of the XY normals are doubled (in the GC and the mean table),
so the PON is effectively female; the PAR is left as is (an XY normal has two PAR copies). XX, ambiguous and unknown
normals are not doubled: one warning per ambiguous or unknown normal, a summary warning when any is ambiguous, and a
warning when main's rule (normalised chrX mean < 75) would double other normals. `CNV/pon.json` records the calls:
`"sex": [{"column": "Cov1", "name": <PON_list.txt entry or null>, "chrx_ratio": 1.0571, "call": "XX", "doubled": false},
…]` and `"sex_thresholds"` (thresholds, PARs, `pon_call` and `doubling` rules). `name` comes from
`<--ponmatrix>/PON_list.txt` when it exists and has one entry per column. `--legacy-stats` keeps main's rule
(normalised chrX mean < 75, all chrX bins) and only records the calls. Check the calls and warnings of a new kit or PON
before using it.

**Samples** (`sample`, `combine`, `sex`): `sex.tsv` from the raw coverage table, i.e. the log2ratio of the normal and
the tumor against the (effectively female) PON, one row per sample column (tumor-only: one row, the tumor's). An empty
table, or one without a finite autosomal log2ratio, fails (exit 2); a table without chrX bins gives `unknown`:

| column | content |
|---|---|
| `sample` | label: `--normal-name`/`--tumor-name` (or from the legacy sample ID; `sex` without names: `normal`/`tumor`; tumor-only: `--tumor-name`, else `--sample`, else `tumor`) |
| `role` | `normal` or `tumor` |
| `column` | `log2ratio1` (normal = sample 1) or `log2ratio2` (tumor); tumor-only: `log2ratio1` (the tumor) |
| `chrx_ratio`, `log2_diff` | `2 ** log2_diff` and median(chrX non-PAR) − median(autosomes), 4 decimals (empty if undefined) |
| `n_chrx_bins`, `n_auto_bins` | finite values used |
| `call` | `XX`, `XY`, `ambiguous` or `unknown` |
| `ascat_gender` | the gender ASCAT uses, the same on every row: the normal's call (a tumor-only run: the tumor's); ambiguous or unknown → `XX`, with a warning |

The pipeline's ASCAT step (`ASCAT_PONCNV`) passes `ascat_gender` to `ascat.loadData(gender = …)`. The sheet sex is
only a cross-check there: a mismatch with `call` is written to `status.tsv`; for an ambiguous or unknown call the
template may use the sheet sex instead of the `XX` fallback (design §3, note 6). A tumor that lost one X is called XY
on its own row, but a pair keeps the normal's gender.

## Failures instead of empty outputs

A task that cannot produce a meaningful result exits non-zero (2 for data problems, 1 for a failing shell pipeline)
with one log line that names the cause. Examples: a missing or mistyped `--gc`/`--genmap`, a track without data for a
chromosome of the sample, an empty pileup, a pileup with another number of samples than the mode needs (2 for a pair,
1 with `--tumor-only`), a BED that is unsorted or does not match the pileup, a sample binned with
another BED or coordinate mode than the CNV PON (checked with `CNV/pon.json` and on the bins themselves), fewer than
90 % of the sample bins found in the PON (`--min-pon-match`), PON chromosomes missing from the sample (a partial
`combine`; `--allow-missing-chroms`), or no SNP left for ASCAT (`--allow-empty` writes header-only tables instead).

## Legacy compatibility

- `--legacy-coords` (`pon-coverage`, `tn2cnv`, `combine`, `sample`) reproduces the old bin labels and ExonPos. Use it
  only to combine new samples with PON coverage files built by the legacy tools; otherwise rebuild the CNV PON from
  `matrix/`. The dev PON `somvar_HAEv7_NovaSeq_v1` has legacy `cov/` and `CNV/` tables built on the **padded**
  SS_HAEv7 BED: run samples against it with `--bed SS_HAEv7_hg38_Padded.bed --legacy-coords`.
- `pon-combine --legacy-stats` reproduces the old PON tables (filter stats not recomputed, all chrX bins doubled by
  main's rule: normalised chrX mean < 75).
- All deliberate changes and the measured deltas on a study sample: `CHANGELOG.md`.

## Layout

`src/mycnv/`: `raw` (tn2cnv/tumor2cnv/pon2cnv), `sites` (tumor-only SNP sites), `bed` (BED reading, bin ExonPos), `tracks` (GC/genmap reading), `trackgen` (GC/genmap building), `combine` (PON and
sample combine), `ponmeta` (PON provenance), `sex` (coverage sex calls), `rolling`, `rolling_cov`, `rolling_snp`, `convert` (ASCAT input),
`pipeline` (steps used by the CLI), `cli`, `plot*`, `cluster` (GMM-caller groundwork; needs scikit-learn), `mawk/`.
Notebooks in `nb/` are the development history and still use the legacy `code/` imports (plan 2.3 moves them to
`import mycnv`). Tests: `tests/README.md`. The subcommands and options are pinned by a snapshot test
(`tests/data/cli_contract.json`, plan C1), because the kumo-variant-calling-nf modules call them with fixed flags. After
a deliberate CLI change, rerun `python tests/make_cli_contract.py` and review the modules that use the changed flags.
