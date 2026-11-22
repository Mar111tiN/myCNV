# mycnv tests

`python -m pytest -q` from the repo root (env: `environment.yml` or the workbench `somvar` env).

## Fixtures

| Path | What | In git |
|---|---|---|
| `data/synth/` | synthetic clean T/N pileups (`chr7.cpile` on HAEv7 covered regions, `multi.cpile` on 3 chromosomes), a synthetic 5-normal PONmatrix, `multi.bed` (`synth.py`) | yes |
| `data/tracks/` | chr7 slices of the public GC (`gc100-10`) and genmap tracks (`make_fixtures.py`) | yes |
| `data/haev7_{covered,padded}_chr7.bed` | SS_HAEv7 BED coordinates (3 columns) over the first 300 kb of chr7 | yes |
| `data/legacy/` | `main` outputs on synthetic inputs (`legacy_outputs.py`) | yes |
| `data/restricted/` | **patient-derived**: chr7 PONmatrix head and cov/snp heads of the somvar HAEv7 PON (28 normals), NHL1 `raw/` slices, the legacy outputs made from them | **no** (git-ignored, plan D14) |

Tests that need `data/restricted/` are skipped when it is absent (`pytest -rs` lists them). To recreate it on the
workbench (dev-space test data synced to `/home/GOFQR/tools/somvar-dev/data`, see `make_fixtures.py`):
```bash
python tests/make_fixtures.py /home/GOFQR/tools/somvar-dev/data          # somvar env
python tests/synth.py                                                     # somvar env
git worktree add --detach /tmp/myCNV-main main
conda run -n somvar-legacy python tests/legacy_outputs.py /tmp/myCNV-main  # legacy env, unmodified main
```

## Parity scripts (not pytest)

- `study_rolling.py legacy|new …`: whole-sample rolling step on an old study sample; `compare_tables.py` compares the
  result with the study's `data/*.cnv.snp.gz`. Results for NHL1_A-B are in `CHANGELOG.md`.
