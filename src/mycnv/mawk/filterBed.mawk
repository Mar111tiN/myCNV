#!/usr/bin/env bash
# filterBed v2.1.0 (myCNV kumodev)
# Filters a position stream (header required, columns 1 = Chr, 2 = Pos) to the positions of a BED file and
# optionally appends a per-chromosome cumulative exonic coordinate (ExonPos).
#
# USAGE:  <stream with header> | bash filterBed.mawk <bedfile> [-x] [-c chrom[,chrom...]] [-L]
#   -x | --use-exon-pos     append ExonPos
#   -c | --chrom            only keep these chromosomes ("7" or "chr7", comma-separated); default: all BED chromosomes
#   -L | --legacy-coords    reproduce the v1.0.1 ExonPos formula (BEDSUM - BEDSTART + pos). It is not
#                           monotonic within a chromosome; use it only together with PON coverage files that were
#                           built by the legacy tools.
#
# Changes vs v1.0.1 (shell/filterBed.mawk on main):
#   - bash shebang (the argument loop is bash syntax; under dash every option was silently ignored)
#   - BED header lines (browser/track/#) are skipped. Each BED chromosome must be one block sorted by start
#     (`sort -k1,1 -k2,2n`; the blocks may come in any order): a chromosome that comes back after another one, or a
#     start below the previous one, stops the tool with an error (v1.0.1 dropped or misassigned positions silently).
#     The stream's chromosomes may come in any order; positions must be sorted within each chromosome.
#   - -c takes a comma-separated list
#   - no early exit after the last BED region: later chromosomes of the stream are no longer dropped
#   - ExonPos = cumulative covered bases before the region + (pos - BEDstart + 1), i.e. contiguous and
#     strictly increasing within a chromosome (v1.0.1 used BEDSUM - BEDSTART + pos)
# Kept on purpose: a BED region [start, end) covers the stream positions start..end (both inclusive), as before,
# so ExonPos keys stay comparable between samples and PONs processed by either version.

useExonicCoords=0
legacyCoords=0
filterChrom=""
PARAMS=""
while (( "$#" )); do
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--use-exon-pos)
        useExonicCoords=1
        shift
        ;;
        -L|--legacy-coords)
        legacyCoords=1
        shift
        ;;
        -c|--chrom)
        if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
            filterChrom=$2
            shift 2
        else
            echo "<filterBed> Error: chromosome argument is missing [-c|--chrom (default=all)]" >&2
            exit 1
        fi
        ;;
        -*|--*=)
        echo "<filterBed> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *)
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done
eval set -- "$PARAMS"
bedFile=$1
if [ -z "$bedFile" ] || [ ! -r "$bedFile" ]; then
    echo "<filterBed> Error: BED file '$bedFile' not readable" >&2
    exit 1
fi

# the BED file is read first (FNR==NR), then the stream from stdin
mawk -v useExonicCoords="$useExonicCoords" -v legacyCoords="$legacyCoords" -v filterChrom="$filterChrom" -v bedFile="$bedFile" '
BEGIN {
    FS = "\t";
    nKeep = (filterChrom != "") ? split(filterChrom, KEEPLIST, ",") : 0;
    for (k = 0; k++ < nKeep;) {
        if (KEEPLIST[k] !~ /^chr/) KEEPLIST[k] = "chr" KEEPLIST[k];
        KEEP[KEEPLIST[k]] = 1;
    }
    printf("<filterBed> Reading bed file %s%s\n", bedFile, nKeep ? " (chromosomes " filterChrom ")" : "") > "/dev/stderr";
    bedCount = 0;
    chromCount = 0;
    lastBedChrom = "";
}
######## BED ########
FNR == NR {
    if ($0 ~ /^(browser|track|#)/ || $0 == "") next;
    # every chromosome one block, sorted by start: otherwise positions would be dropped or misassigned silently
    if ($1 != lastBedChrom) {
        if ($1 in BEDSEEN) {
            printf("<filterBed> Error: %s line %s: %s comes back after %s; sort the BED (sort -k1,1 -k2,2n)\n", bedFile, FNR, $1, lastBedChrom) > "/dev/stderr";
            exit 1;
        }
        BEDSEEN[$1] = 1;
        lastBedChrom = $1;
    } else if ($2 + 0 < lastBedStart) {
        printf("<filterBed> Error: %s line %s: %s:%s is listed after start %s; sort the BED (sort -k1,1 -k2,2n)\n", bedFile, FNR, $1, $2, lastBedStart) > "/dev/stderr";
        exit 1;
    }
    lastBedStart = $2 + 0;
    if (nKeep && !($1 in KEEP)) next;
    bedCount++;
    if (!($1 in CHROMSTART)) {
        # first region of a new chromosome
        if (chromCount) CHROMEND[CHROMNAME[chromCount]] = bedCount - 1;
        CHROMNAME[++chromCount] = $1;
        CHROMSTART[$1] = bedCount;
        exonicCoord = 0;
    }
    BEDSTART[bedCount] = $2 + 0;
    BEDEND[bedCount] = $3 + 0;
    exonicCoord += $3 - $2 + 1;
    BEDSUM[bedCount] = exonicCoord;
    next;
}
######## STREAM HEADER ########
!headerDone {
    if (chromCount) CHROMEND[CHROMNAME[chromCount]] = bedCount;
    headerDone = 1;
    currentChrom = "";
    if ($0 ~ /^Chr\t/) {
        printf("%s", $0);
        if (useExonicCoords) printf("\tExonPos");
        printf("\n");
        next;
    }
    print("<filterBed> No header detected") > "/dev/stderr";
}
######## DATA ########
{
    chrom = $1;
    if (chrom != currentChrom) {
        currentChrom = chrom;
        if (!(chrom in CHROMSTART)) { skipChrom = 1; next; }
        skipChrom = 0;
        bedPointer = CHROMSTART[chrom];
        lastPointer = CHROMEND[chrom];
    }
    if (skipChrom) next;
    pos = $2 + 0;
    while (pos >= BEDSTART[bedPointer]) {
        if (pos > BEDEND[bedPointer]) {
            if (++bedPointer > lastPointer) { skipChrom = 1; next; }
        } else {
            printf("%s", $0);
            if (useExonicCoords) {
                if (legacyCoords) {
                    exonicCoord = BEDSUM[bedPointer] - BEDSTART[bedPointer] + pos;
                } else {
                    exonicCoord = BEDSUM[bedPointer] - BEDEND[bedPointer] + pos;
                }
                printf("\t%s", exonicCoord);
            }
            printf("\n");
            break;
        }
    }
}' "$bedFile" -
