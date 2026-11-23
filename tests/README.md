# mycnv tests

`python -m pytest -q` from the repo root (env: `environment.yml` or the workbench `somvar` env).

## Tiers

- **Public** (always runs, also in CI and `docker build --target test`): mawk tools, bins and coordinates
  (`test_coords.py`), tracks and failure modes (`test_tracks.py`, `test_cli.py`), and value parity with `main` for
  every step on synthetic inputs (`test_parity_synth.py`, `test_mawk.py`, `test_raw.py`), the coverage sex calls, the
  PON chrX doubling and `sex.tsv` (`test_sex.py`), the CLI end to end, the packaging contract (`test_packaging.py`).
- **Fixture smoke** (`tests/smoke/smoke_mycnv.sh tests/data/plumbing`): the CLI end to end on the non-secret tier-A
  fixture (`tests/data/plumbing/README.md`); CI runs it inside the runtime image. Not collected by pytest.
- **Restricted** (skipped without `data/restricted/`; `pytest -rs` lists them): the same parity checks on slices of the
  real HAEv7 PON and of the NHL1 study sample (`test_pon.py`, `test_combine.py`, `test_rolling.py`, `test_convert.py`,
  parts of `test_mawk.py` / `test_raw.py`).

## Fixtures

| Path | What | In git |
|---|---|---|
| `data/synth/` | synthetic clean T/N pileups (`chr7.cpile` on HAEv7 covered regions, `multi.cpile` on 3 chromosomes), PONmatrix files (`matrix/`, `matrix_multi/`), `multi.bed`, its GC/genmap tracks (`tracks/`), a raw mpileup snippet, a synthetic CNV PON (`pon/`, 8 normals, 2 male; chrX at 3.0 Mb, outside PAR1) and T/N sample (`sample/`) (`synth.py`) | yes |
| `data/tracks/` | chr7 slices of the public GC (`gc100-10`) and genmap tracks (`make_fixtures.py`) | yes |
| `data/haev7_{covered,padded}_chr7.bed` | SS_HAEv7 BED coordinates (3 columns) over the first 300 kb of chr7 | yes |
| `data/legacy/` | `main` outputs on the synthetic inputs (`legacy_outputs.py`) | yes |
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

## Parity scripts (not pytest)

- `study_rolling.py legacy|new …`: whole-sample rolling step on an old study sample; `compare_tables.py` compares the
  result with the study's `data/*.cnv.snp.gz`. Results for NHL1_A-B are in `CHANGELOG.md`.
