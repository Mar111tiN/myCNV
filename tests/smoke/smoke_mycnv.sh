#!/usr/bin/env bash
# Smoke test of the mycnv runtime image on fixture tier A (run INSIDE the container).
#   podman run --rm --platform linux/amd64 -v "$PWD/tests/data/plumbing:/data:ro" -v "$PWD/tests/smoke:/smoke:ro" \
#       localhost/mycnv:<version> bash /smoke/smoke_mycnv.sh /data
# Inputs are the ebscore outputs in /data/expected/ebscore (PONmatrix/matrix + pair.pileup.gz), so this test does not
# need the ebscore image. Writes to /tmp/smoke inside the container; compares with /data/expected/mycnv.
# The JPGs are only checked for existence. Exit 0 = all PASS.
set -uo pipefail

D=${1:-/data}
S=$(cd "$(dirname "$0")" && pwd)
I=$D/expected/ebscore
E=$D/expected/mycnv
O=/tmp/smoke/mycnv
P=fixture_tumor-normal
mkdir -p "$O"
fail=0
step() { echo "== $*"; }
check() { if "$@"; then return 0; else fail=1; return 1; fi; }
same() { check cmp "$1" "$2" && echo "PASS $(basename "$2") identical"; }
table() { check python "$S/compare_tables.py" "$@"; }

step "environment"
echo "user=$(id -un) uid=$(id -u)  (must not be root)"
[ "$(id -u)" != "0" ] || { echo "FAIL runs as root"; fail=1; }
mycnv --version
samtools --version | head -1
mawk -W version 2>&1 | head -1
# `mycnv tracks genmap` needs genmap (environment.yml); no `check` in a pipeline: it would run in a subshell
if genmap --version >/dev/null 2>&1; then genmap --version | head -1; else echo "FAIL genmap missing"; fail=1; fi
python -c "import pandas, numpy, matplotlib; print('pandas', pandas.__version__, 'numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)"
echo "MPLBACKEND=${MPLBACKEND:-unset} MPLCONFIGDIR=${MPLCONFIGDIR:-unset}"
ps -o rss= -p $$ >/dev/null && echo "ps OK" || { echo "FAIL ps (procps) missing"; fail=1; }

step "pon-coverage"
check mycnv pon-coverage --ponmatrix "$I/PONmatrix" --bed "$D/bed/testkit.padded.bed" \
    --gc "$D/static/gc100-10" --genmap "$D/static/genmap" -o "$O/PONmatrix" --quiet
table "$E/PONmatrix/cov/chr21.cov.gz" "$O/PONmatrix/cov/chr21.cov.gz"
table "$E/PONmatrix/snp/chr21.snp.gz" "$O/PONmatrix/snp/chr21.snp.gz"
same "$E/PONmatrix/cov/chr21.cov.json" "$O/PONmatrix/cov/chr21.cov.json"

step "pon-combine"
check mycnv pon-combine --ponmatrix "$O/PONmatrix" -o "$O/PONmatrix" --chrom chr21 --quiet
for f in pon.cov.full pon.cov.fullGC pon.cov.filter pon.cov.filterGC pon.snp; do
    table "$E/PONmatrix/CNV/$f.gz" "$O/PONmatrix/CNV/$f.gz"
done
same "$E/PONmatrix/CNV/pon.json" "$O/PONmatrix/CNV/pon.json"

step "sample (tumor/normal pair, with plots)"
check mycnv sample --pileup "$I/pair.pileup.gz" --ponmatrix "$O/PONmatrix" --bed "$D/bed/testkit.padded.bed" \
    --gc "$D/static/gc100-10" --genmap "$D/static/genmap" --sample $P --tumor-name tumor --normal-name normal \
    -o "$O/pair" --plots --allow-empty --quiet
table "$E/pair/raw/$P.cov.gz" "$O/pair/raw/$P.cov.gz"
table "$E/pair/raw/$P.snp" "$O/pair/raw/$P.snp"
table "$E/pair/data/$P.cnv.gz" "$O/pair/data/$P.cnv.gz"
table "$E/pair/data/$P.cnv.snp.gz" "$O/pair/data/$P.cnv.snp.gz"
# coverage sex calls (M5): the chr21-only fixture has no chrX bins, so both rows are unknown and ASCAT gets XX
table "$E/pair/sex.tsv" "$O/pair/sex.tsv"
for f in "$E"/pair/pre/*.tsv; do
    table "$f" "$O/pair/pre/$(basename "$f")"
done
# the JPGs are not stored in expected/ (only their existence is checked): the 6 plots of a pair whose SNP table is empty
for name in normal.cnv normal.cov.raw normal.snp.raw tumor.cnv tumor.cov.raw tumor.snp.raw; do
    f="$O/pair/plots/${P}_${name}.jpg"
    if [ -s "$f" ]; then
        echo "PASS $(basename "$f") written"
    else
        echo "FAIL $(basename "$f") missing or empty"; fail=1
    fi
done

step "sample without --allow-empty must exit 2 (no SNP passes the filter on this fixture)"
mycnv sample --pileup "$I/pair.pileup.gz" --ponmatrix "$O/PONmatrix" --bed "$D/bed/testkit.padded.bed" \
    --gc "$D/static/gc100-10" --genmap "$D/static/genmap" --sample $P --tumor-name tumor --normal-name normal \
    -o "$O/pair_strict" --quiet 2>/dev/null
rc=$?
if [ "$rc" = 2 ]; then echo "PASS exit 2"; else echo "FAIL exit $rc, expected 2"; fail=1; fi

step "tracks gc (samtools faidx inside the image; must reproduce the fixture GC slice)"
check mycnv tracks gc --fasta "$D/ref/chr21.fa.gz" --chrom chr21 -o "$O/tracks/gc100-10" --quiet
check python "$S/compare_track_slice.py" "$D/static/gc100-10/chr21.gc100-10.gz" "$O/tracks/gc100-10/chr21.gc100-10.gz"

step "error handling (missing input must exit non-zero)"
if mycnv pon-combine --ponmatrix /nonexistent -o "$O/x" --quiet 2>/dev/null; then
    echo "FAIL missing input exited 0"; fail=1
else
    echo "PASS missing input exits non-zero"
fi

echo
if [ "$fail" = 0 ]; then echo "ALL PASS (mycnv)"; else echo "SOME CHECKS FAILED (mycnv)"; fi
exit $fail
