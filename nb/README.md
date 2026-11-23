# myCNV notebooks

Development history of the CNV code. On `kumodev` the legacy `code/` modules live in the `mycnv` package
(`src/mycnv/`), so a notebook only runs if it was ported to `import mycnv` (plan 2.3, branch `dev/nb-port`).
`main` keeps every notebook unchanged with its outputs.

## Status

| Notebook | Status | Inputs (via `config/config_dev.yaml`) | Writes (`<out_root>/<sample>/`) |
|---|---|---|---|
| `0A_getPONdata` | ported | PON `matrix/{chrom}.pon.gz`, BED, GC, genmap (chromosomes: `pon_build.chroms`) | `<out_root>/PON/{cov,snp}/{chrom}.*.gz` |
| `0B_combinePONdata` | ported | PON `cov/`, `snp/` (chromosomes: `pon_build.combine_chroms`) | `<out_root>/PON/CNV/pon.cov.{full,fullGC,filter,filterGC}.gz`, `pon.snp.gz`, `PON/figures/` |
| `1_getCOV+SNP` | not ported | BAMs or paired pileups (none exist for the study data) | — |
| `2_combineSampledata` | not ported | the per-chromosome output of notebook 1 | — (its result is the study's `raw/`) |
| `3_rollingCoverage` | ported | study `raw/{pair}.cov.gz` | `cov/{pair}.roll.cov.gz`, `plot/{pair}.log2gauss.jpg` |
| `4_rollingSNP` | ported | study `raw/{pair}.snp` (plain text) + notebook 3's `roll.cov.gz` | `data/{pair}.cnv.snp.gz`, `data/{pair}.cnv.gz`, `figures/` |
| `5_CNV2ASCAT` | ported | `data/{pair}.cnv.snp.gz` (`chain_from`: this run or the study) | `pre/{pair}_{N,T,germline}_{baf,logr}.tsv` |
| `visualize` | ported (cluster section skipped: Gen-1 `*.dbscan.csv` missing) | `data/{pair}.cnv.gz` and `.cnv.snp.gz` | nothing (plots only) |
| `6_clusterCNV`, `AML_*`, `CNVinspect`, `GMMcenter`, `PMBL_rollingCoverage`, `visualizeCluster` | not ported (Gen-1) | the deleted `codeCNV` package and Gen-1 intermediates (`.cluster`, `.roll.snp`, `.dbscan.csv`) | — |
| `gisticLike` | not ported | external `~/Sites/Bio/{pytools,pycnv}` and `CNV_coll3.txt` | — |

Verified 2026-11-23 (headless, kernel `somvar`, mycnv `kumodev` `febd0c3`): `0A`/`0B` on the D21 HAEv7 PON, `3 → 4 → 5 →
visualize` on NHL1_A-B and 03_A-B run without a cell error.
- `0A` chr21 equals the CLI build (`mycnv pon-coverage`) in every shared column; the notebook adds the genmap columns.
- `0B` on 3 chromosomes differs from the genome-wide CLI build by one scale factor per normal (normalisation over fewer
  chromosomes), as expected.
- NHL1: the same 12,939 SNPs as the study; AML 03: 51 more SNPs, explained in `../CHANGELOG.md` (fallSNP max-scaling on
  the legacy coordinates). Unmodified `main` reproduces both study files exactly.

What the port changed in the ported notebooks:
- `sys.path.append('../code')` + module imports → `import mycnv` (`mycnv.rolling`, `rolling_cov`, `rolling_snp`,
  `combine`, `convert`, `raw`, `plot`, `plot_helpers`); renamed functions: `mergeCov2SNP` → `merge_cov_to_snp`,
  `get_fallSNP` → `get_fall_snp`, `remergeCNV` → `remerge_cnv`, `write_ASCAT` → `write_ascat`,
  `combine_PON_coverage`/`make_PON_coverage`/`make_PON_snp` → `combine_pon_coverage`/`make_pon_coverage` (returns a
  dict)/`make_pon_snp`, `normalize_GC` → `normalize_gc`, `PON2CNV` → `pon2cnv`.
- Inline copies of module functions (the "step by step" cells) are replaced by imports of the canonical `mycnv`
  function; the cells that call them are unchanged. The behaviour changes of `mycnv` vs `main` are in `../CHANGELOG.md`.
- Hard-coded Dropbox paths and the sample `03_A-B` → `dev_config.load_dev()` (`../config/config_dev.yaml`).
- Explicit `import numpy as np, pandas as pd` (main relied on an IPython startup file). `pd.set_option(
  'mode.chained_assignment', None)` is dropped (no effect with pandas 3).
- Cells that only worked in the author's interactive session got a minimal fix, marked with a `# main …` comment:
  `3`: `df.drop(['L', 'R'])` of columns that never exist on a fresh run; `0B`: `normGC_cov_df` used before it was
  defined, undefined `cov_agg`, `amazonize`/`compute_stats`/`remove_outliers` returned lists.
- `visualize`: the Gen-1 `plot_snp` with coverage overlay is renamed `plot_snp_cluster` (it shadowed `plot_snp`).
- `3`: the plot dicts `llh`/`llhdiff` are `llh_plot`/`llhdiff_plot`; the step-by-step cells rebind `llh` to the
  function, so main's region plot failed on a fresh run (`TypeError: 'function' object is not subscriptable`).
- Outputs are cleared in git. Never commit executed copies: the study data are consented, workbench-only data.

## Running

Install `mycnv` into the kernel's env (`pip install -e ..`; the workbench `somvar` env has the `kumodev` checkout
installed editable, so set `PYTHONPATH=<this worktree>/src` to run another branch). Then, from this folder:

```bash
# interactive: kernel "Python (somvar)"; the sample is `sample:` in config/config_dev.yaml
# headless, one pair (keys under `samples:`), outputs to <out_root>/<key>/
MYCNV_NB_SAMPLE=D1 jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=somvar \
    --ExecutePreprocessor.timeout=3600 --output-dir <out_root>/D1 3_rollingCoverage.ipynb
```
Order per pair: 3 → 4 → 5 → visualize. Overrides: `MYCNV_NB_CONFIG` (another config file), `MYCNV_NB_DATA` (data
root), `MYCNV_NB_OUT` (output root). Notebooks 3 and 4 hold the whole sample in memory (~9 M coverage bins); plan
on ~8 GB per kernel.
