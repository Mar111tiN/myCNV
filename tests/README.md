# mycnv tests

`python -m pytest -q` from the repo root (env: `environment.yml` or the workbench `somvar` env).

## Tiers

- **Public** (always runs, also in CI and `docker build --target test`): mawk tools, bins and coordinates
  (`test_coords.py`), tracks and failure modes (`test_tracks.py`, `test_cli.py`), and value parity with `main` for
  every step on synthetic inputs (`test_parity_synth.py`, `test_mawk.py`, `test_raw.py`), the coverage sex calls, the
  PON chrX doubling and `sex.tsv` (`test_sex.py`), tumor-only samples (`test_tumor_only.py`: het SNPs at the sites
  within the tumor VAF band and depth, the sites file rules (formats, headers, multi-allelic sites), one-sample tables
  equal to a pair's tumor columns, mock germline, one-row `sex.tsv`, pileups on a pipe, CLI and error cases), the CLI
  end to end, the packaging contract (`test_packaging.py`), the CLI contract (`test_cli_contract.py`, below).
- **Fixture smoke** (`tests/smoke/smoke_mycnv.sh tests/data/plumbing`): the CLI end to end on the non-secret tier-A
  fixture (`tests/data/plumbing/README.md`); CI runs it inside the runtime image. Not collected by pytest.
- **Restricted** (skipped without `data/restricted/`; `pytest -rs` lists them): the same parity checks on slices of the
  real HAEv7 PON and of the NHL1 study sample (`test_pon.py`, `test_combine.py`, `test_rolling.py`, `test_convert.py`,
  parts of `test_mawk.py` / `test_raw.py`).

## Fixtures

| Path | What | In git |
|---|---|---|
| `data/synth/` | synthetic clean T/N pileups (`chr7.cpile` on HAEv7 covered regions, `multi.cpile` on 3 chromosomes), PONmatrix files (`matrix/`, `matrix_multi/`), `multi.bed`, its GC/genmap tracks (`tracks/`), a raw mpileup snippet, a synthetic CNV PON (`pon/`, 8 normals, 2 male; chrX at 3.0 Mb, outside PAR1) and T/N sample (`sample/`) (`synth.py`) | yes |
| (test time) | tumor-only inputs made from the files above in `tmp_path`: the one-sample pileup = the tumor columns of `multi.cpile` (`synth.make_tumor_pileup`), a sites VCF/TSV with its het, homozygous, reference-like and absent positions (`synth.make_snp_sites`, `synth.sites_vcf`), one-sample tables = the tumor columns of `sample/` | no (made at test time) |
| `data/tracks/` | chr7 slices of the public GC (`gc100-10`) and genmap tracks (`make_fixtures.py`) | yes |
| `data/haev7_{covered,padded}_chr7.bed` | SS_HAEv7 BED coordinates (3 columns) over the first 300 kb of chr7 | yes |
| `data/legacy/` | `main` outputs on the synthetic inputs (`legacy_outputs.py`) | yes |
| `data/cli_contract.json` | snapshot of the CLI's subcommands and options (`make_cli_contract.py`, plan C1) | yes |
| `data/restricted/` | **patient-derived**: chr7 PONmatrix head and cov/snp heads of the somvar HAEv7 PON (28 normals), NHL1 `raw/` slices, the legacy outputs made from them | **no** (git-ignored, plan D14; not in the Docker build context) |

Regenerate (workbench; dev-space test data synced to `/home/GOFQR/tools/somvar-dev/data`, see `make_fixtures.py`;
scratch outside `/tmp`):
```bash
python tests/make_fixtures.py /home/GOFQR/tools/somvar-dev/data          # somvar env (restricted + tracks)
python tests/synth.py                                                     # somvar env (synthetic; reproducible bytes)
git worktree add --detach <scratch>/myCNV-main main
conda run -n somvar-legacy python tests/legacy_outputs.py <scratch>/myCNV-main  # legacy env, unmodified main
```
`legacy_outputs.py` writes the public outputs always and the restricted ones only if `data/restricted/` exists.

## CLI contract (plan C1)

The kumo-variant-calling-nf modules call `mycnv` with fixed flags (design §3 process table: `MYCNV_PONCOVERAGE` →
`pon-coverage`, `MYCNV_PONCOMBINE` → `pon-combine`, `MYCNV_SAMPLE` → `sample`, paired and `--tumor-only`), so
`test_cli_contract.py` compares `mycnv.cli.build_parser()` with the snapshot `data/cli_contract.json` (schema
`mycnv-cli-contract/1`, the same layout as ebscore's). For every subcommand, the nested `tracks gc` and `tracks genmap`
included, it holds the sorted option strings, and per option the flags, required, the action kind (`store`,
`store_true`, `append`, `extend`, …), `nargs`, `type`, `choices` and the default when it is a constant; also the
mutually exclusive groups (`pon-coverage`: `--ponmatrix` | `--matrix`, required), positionals, and the `--version`
output format with the version masked (`mycnv <version>`; versions.yml parses it). Help texts and metavars are not part
of it (argparse formats `--help` differently between Python versions; the snapshot is the same on Python 3.9 and 3.15).
The comparison is type-exact (a default `1` → `1.0` or `0` → `false` is a change). The test also parses the design §3
module command lines with prefix matching off and runs the CLI's cross-option checks on them (e.g. `--tumor-only`
needs `--snp-sites`), and checks that no (sub)command parser allows prefix matching.

On a mismatch the message lists each changed entry and the modules of the changed subcommands. After a deliberate CLI
change:
```bash
python tests/make_cli_contract.py   # rewrites the JSON from this checkout's src/, prints the changes and the modules
```
then review the kumo-variant-calling-nf modules that call the changed flags (design §3) and commit the JSON with the
change. The message names the `mycnv` the test imported. If it warns about a stale install (another copy than this
checkout's `src/`, with other code: on the workbench the `somvar` env's editable install points at the `kumodev`
checkout, so a worktree needs `PYTHONPATH=src`), do not regenerate: set `PYTHONPATH=src` or reinstall.

## Parity scripts (not pytest)

- `study_rolling.py legacy|new …`: whole-sample rolling step on an old study sample; `compare_tables.py` compares the
  result with the study's `data/*.cnv.snp.gz`. Results for NHL1_A-B are in `CHANGELOG.md`.
