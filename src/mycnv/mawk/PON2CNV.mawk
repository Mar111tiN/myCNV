#!/usr/bin/env bash
# PON2CNV v2.0.0 (myCNV kumodev)
# Computes binned per-sample coverage and pooled SNP VAFs from a PONmatrix (ebscore pile2count format).
#
# INPUT (optionally | filterBed -x):
#   Chr  Start  Ref  A  G  C  T  I  D  Depth  [ExonPos]      cells: s1|s2|..|sN=s1|s2|..|sN (+ strand = - strand)
# USAGE:
#   gunzip < matrix/chr7.pon.gz | bash filterBed.mawk covered.bed -x -c chr7 | bash PON2CNV.mawk -x -o snp/chr7.snp
# OPTIONS
#   -c | --min-coverage     <INT=0>      minimum (binned) coverage of any sample for coverage output
#   -x | --keep-exon-pos                 keep the ExonPos column (if present)
#   -b | --bin-size         <INT=10>     bin size for the mean coverage
#   -o | --output-snp-file  <path>       SNP output file (required; ".snp" is appended if missing)
#   -v | --min-vaf          <FLOAT=0.25> minimum pooled VAF for SNP output (exclusive)
#   -d | --min-depth        <INT=50>     minimum pooled depth for SNP output
#   -L | --legacy-bins                   reproduce the v1.0.0 bin labels (see pile2CNV.mawk)
# OUTPUT
#   coverage (stdout): Chr Pos [ExonPos] Cov1 .. CovN
#   SNPs (-o file):    Chr Pos [ExonPos] VAF Depth            (pooled over all PON samples and both strands)
#
# Changes vs v1.0.0 (shell/PON2CNV.mawk on main):
#   - bash shebang; per-chromosome binning and bin labels as in pile2CNV v2.0.0
#   - the first base of every bin after the first is counted. v1.0.0 reset the bin sum with the undefined
#     $COL[s] (= $0, numeric 0), so that base was lost (about 1/binSize of the coverage).

keepExonPos=0
legacyBins=0
PARAMS=""
while (( "$#" )); do
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--keep-exon-pos) keepExonPos=1; shift ;;
        -L|--legacy-bins) legacyBins=1; shift ;;
        -o|--output-snp-file|-c|--min-coverage|-b|--bin-size|-v|--min-vaf|--min-VAF|-d|--min-depth)
        if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
            case "$1" in
                -o|--output-snp-file) snpFile=$2 ;;
                -c|--min-coverage) minCov=$2 ;;
                -b|--bin-size) binSize=$2 ;;
                -v|--min-vaf|--min-VAF) minVAF=$2 ;;
                -d|--min-depth) minDepth=$2 ;;
            esac
            shift 2
        else
            echo "<PON2CNV> Error: argument for $1 is missing" >&2
            exit 1
        fi
        ;;
        -*|--*=) echo "<PON2CNV> Error: Unsupported flag $1" >&2; exit 1 ;;
        *) PARAMS="$PARAMS $1"; shift ;;
    esac
done
if [ -z "$snpFile" ]; then
    echo "<PON2CNV> Error: SNP output file is required [-o|--output-snp-file]" >&2
    exit 1
fi

mawk -v keepExonPosArg="$keepExonPos" -v legacyBins="$legacyBins" -v minCov="${minCov-0}" -v binSize="${binSize-10}" \
     -v minVAF="${minVAF-0.25}" -v minDepth="${minDepth-50}" -v snpFile="$snpFile" '
BEGIN { FS = "\t"; }
NR == 1 {  ### HEADER
    hasExonPos = ($NF == "ExonPos");
    keepExonPos = (hasExonPos && keepExonPosArg);
    minCov += 0; binSize += 0; minVAF += 0; minDepth += 0;
    if (minCov) covThresh = minCov * binSize - 1;
    printf("<PON2CNV> Bin coverage: binSize=%s | minCoverage=%s%s\n", binSize, minCov, legacyBins ? " | legacy bin labels" : "") > "/dev/stderr";
    VAFPRECISION = 4;
    PREC = 10 ^ VAFPRECISION;
    printf("<PON2CNV> heteroSNP: minDepth=%s | totalVAF > %s\n", minDepth, minVAF) > "/dev/stderr";
    if (snpFile !~ /\.snp/) snpFile = snpFile ".snp";
    printf("<PON2CNV> Writing snpData to %s\n", snpFile) > "/dev/stderr";
    baseHeader = (keepExonPos) ? "Chr\tPos\tExonPos" : "Chr\tPos";
    printf("%s", baseHeader);
    printf("%s\tVAF\tDepth\n", baseHeader) > snpFile;
    for (f = 0; f++ < NF;) if ($f == "Depth") DEPTHCOL = f;
    haveBin = 0;
    next;
}
NR == 2 {  ### sample count from the first data line, then the coverage header
    split($DEPTHCOL, SPLIT, "=");
    samples = gsub(/\|/, "", SPLIT[1]) + 1;
    print("<PON2CNV> ", samples, "samples detected") > "/dev/stderr";
    zero = "0";
    for (s = 0; s++ < samples;) {
        printf("\tCov%s", s);
        if (s > 1) zero = zero "|0";
        POS[s] = 0;
        NEG[s] = 0;
    }
    zero = zero "=" zero;
    printf("\n");
}

function aggString(col,   PARTS, sum, s) {
    # sum of one cell over samples and strands; leaves the per-sample values in POS/NEG
    sum = 0;
    split(col, PARTS, "=");
    split(PARTS[1], POS, "|");
    split(PARTS[2], NEG, "|");
    for (s = 0; s++ < samples;) sum += POS[s] + NEG[s];
    return sum;
}
function flushBin(label,   s, noOutput) {
    if (minCov) {
        noOutput = 1;
        for (s = 0; s++ < samples;) if (COVSUM[s] > covThresh) { noOutput = 0; break; }
        if (noOutput) return;
    }
    printf("%s\t%s", binChr, label);
    if (keepExonPos) printf("\t%s", binExon);
    for (s = 0; s++ < samples;) printf("\t%s", COVSUM[s] / binSize);
    printf("\n");
}
function startBin(   s) {
    # POS/NEG hold the depths of the current line
    binChr = $1;
    bin = thisBin;
    # legacy: ExonPos extrapolated to the bin start; fixed: ExonPos of the first covered base (monotonic)
    if (keepExonPos) binExon = legacyBins ? $NF - (pos % binSize) : $NF;
    for (s = 0; s++ < samples;) COVSUM[s] = POS[s] + NEG[s];
    haveBin = 1;
}

$DEPTHCOL != zero {  # skip positions without coverage in any sample
    pos = $2;
    ##### pooled VAF
    totalAlt = 0;
    for (col = 3; ++col < DEPTHCOL;) {
        if ($col == zero) continue;
        totalAlt += aggString($col);
    }
    # the Depth cell is aggregated last, so POS/NEG hold the depths below
    totalDepth = aggString($DEPTHCOL);
    totalVAF = int((totalAlt / totalDepth) * PREC) / PREC;
    if ((totalVAF > minVAF) && (totalDepth >= minDepth)) {
        outstring = $1 "\t" pos;
        if (keepExonPos) outstring = outstring "\t" $NF;
        outstring = outstring "\t" totalVAF "\t" totalDepth;
        print(outstring) > snpFile;
    }
    ##### coverage
    thisBin = int(pos / binSize);
    if (!haveBin) {
        startBin();
    } else if ($1 != binChr) {
        flushBin(bin * binSize);
        startBin();
    } else if (thisBin != bin) {
        flushBin(legacyBins ? thisBin * binSize : bin * binSize);
        startBin();
    } else {
        for (s = 0; s++ < samples;) COVSUM[s] += POS[s] + NEG[s];
    }
}
END {
    if (haveBin) flushBin(bin * binSize);
}'
