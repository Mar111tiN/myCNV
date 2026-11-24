# Tier-A plumbing fixture (mycnv copy)

Small end-to-end test data for `tests/smoke/smoke_mycnv.sh`, which CI runs inside the runtime image (as the image user
and as an arbitrary uid). Non-secret; no study or patient data.

| Path | What |
|---|---|
| `ref/chr21.fa.gz` (+ `.fai`, `.gzi`, `.dict`) | GRCh38 chr21 from nf-core/test-datasets, **every base outside chr21:46,099,000-46,301,000 replaced by N**: the test window chr21:46.1-46.3 Mb plus 1 kb on each side (BAQ and the +-50 bp GC windows). Length and coordinates are those of chr21, so the BAM headers and the `.fai` stay valid. The smokes give the same results as with the unmasked chr21 |
| `bed/testkit.padded.bed`, `bed/testkit.window.bed` | the Agilent SureSelect Human All Exon V7 (SS_HAEv7) padded design clipped to the window; the window |
| `static/` | slices of generic hg38 tracks: GC (`gc100-10`, computed from the reference) and genmap mappability |
| `expected/ebscore/` | ebscore outputs on this data (PONmatrix, paired pileup; the inputs of the mycnv smoke) |
| `expected/mycnv/` | mycnv outputs (CNV PON, pair raw/data/pre); the plot JPGs are not stored (the smoke checks only that they are written) |
| `manifest.json` | parameters, sources (public origins), md5 of every file |

Sources and licences: the reads and the chr21 reference come from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets) (MIT licence); the BED slice from the Agilent
SS_HAEv7 design files (publication cleared by the repo owner, 2026-11-23).

Run the smoke without Docker from the repo root (needs the package and its tools on PATH):
`bash tests/smoke/smoke_mycnv.sh tests/data/plumbing`. When an output changes on purpose, regenerate `expected/` in
the same commit.
