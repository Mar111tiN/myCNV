#!/usr/bin/env bash
# pile2CNV v2.1.0 (myCNV kumodev)
# Computes binned coverage and het-SNP data from a clean T/N pileup with depth columns.
#
# INPUT (cleanpileup -d, optionally | filterBed -x), first sample = normal:
#   Chr  Start  Ref  Depth1  Read1  Depth2  Read2 ... [ExonPos]
# USAGE:
#   samtools mpileup -f ref.fa -q 20 -Q 25 -l kit.padded.bed NORMAL TUMOR | bash cleanpileup.mawk -d |
#   bash filterBed.mawk kit.padded.bed -x | bash pile2CNV.mawk -x -o out.snp [options] > out.cov
# OPTIONS
#   -c | --min-coverage     <INT=0>     minimum (binned) coverage of any sample for coverage output
#   -x | --keep-exon-pos                keep the ExonPos column (if present)
#   -b | --bin-size         <INT=10>    bin size for the mean coverage
#   -o | --output-snp-file  <path>      het-SNP output file (required; ".snp" is appended if missing)
#   -v | --min-vaf          <FLOAT=0.1> minimum normal VAF for het-SNP output
#   -V | --max-vaf          <FLOAT=0.9> maximum normal VAF for het-SNP output
#   -d | --min-depth        <INT=10>    minimum depth (every sample) for het-SNP output
#   -L | --legacy-bins                  reproduce the v1.0.1 bin labels (see below)
# OUTPUT
#   coverage (stdout):  Chr Pos [ExonPos] Cov1 Cov2 ...     (mean depth per bin; Pos = bin start)
#   het SNPs (-o file): Chr Pos [ExonPos] Depth1 VAF1 Depth2 VAF2 ...
#
# Changes vs v1.0.1 (shell/pile2CNV.mawk on main):
#   - bash shebang (the argument loop is bash syntax; under dash every option was silently ignored)
#   - streams with several chromosomes are binned per chromosome: v1.0.1 printed the last bin of a chromosome
#     under the next chromosome's name and could merge bins of two chromosomes with the same bin index
#   - Pos is the start of the bin that holds the coverage. v1.0.1 labelled every bin with the start of the next
#     bin it met (the last bin of a chromosome with its own start, which could duplicate the label before it).
#     The bin ExonPos is the ExonPos of the first base of the bin that the stream holds (v1.0.1: that value minus
#     pos % binSize, which is not monotonic at region starts). That value depends on which bases the stream holds,
#     so mycnv (raw.tn2cnv / raw.pon2cnv) replaces it with the ExonPos of the bin's first BED base
#     (mycnv.bed.bin_exon_pos), which samples and PONs share. -L restores both v1.0.1 labels, e.g. to merge with PON
#     coverage files built by the legacy tools.
#   - v2.1.0: het-SNP alt reads include reverse-strand insertions (cleanpileup "i"); v1.0.1 counted A/C/G/T, D/d and
#     only the forward-strand "I", so het insertions got about half their VAF (in both bin modes)

keepExonPos=0
legacyBins=0
PARAMS=""
while (( "$#" )); do
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--keep-exon-pos) keepExonPos=1; shift ;;
        -L|--legacy-bins) legacyBins=1; shift ;;
        -o|--output-snp-file|-c|--min-coverage|-b|--bin-size|-v|--min-vaf|--min-VAF|-V|--max-vaf|--max-VAF|-d|--min-depth)
        if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
            case "$1" in
                -o|--output-snp-file) snpFile=$2 ;;
                -c|--min-coverage) minCov=$2 ;;
                -b|--bin-size) binSize=$2 ;;
                -v|--min-vaf|--min-VAF) minVAF=$2 ;;
                -V|--max-vaf|--max-VAF) maxVAF=$2 ;;
                -d|--min-depth) minDepth=$2 ;;
            esac
            shift 2
        else
            echo "<pile2CNV> Error: argument for $1 is missing" >&2
            exit 1
        fi
        ;;
        -*|--*=) echo "<pile2CNV> Error: Unsupported flag $1" >&2; exit 1 ;;
        *) PARAMS="$PARAMS $1"; shift ;;
    esac
done
if [ -z "$snpFile" ]; then
    echo "<pile2CNV> Error: het-SNP output file is required [-o|--output-snp-file]" >&2
    exit 1
fi

mawk -v keepExonPosArg="$keepExonPos" -v legacyBins="$legacyBins" -v minCov="${minCov-0}" -v binSize="${binSize-10}" \
     -v minVAF="${minVAF-0.1}" -v maxVAF="${maxVAF-0.9}" -v minDepth="${minDepth-10}" -v snpFile="$snpFile" '
BEGIN { FS = "\t"; }
NR == 1 {  ### HEADER
    hasExonPos = ($NF == "ExonPos");
    keepExonPos = (hasExonPos && keepExonPosArg);
    minCov += 0; binSize += 0; minVAF += 0; maxVAF += 0; minDepth += 0;
    if (minCov) covThresh = minCov * binSize - 1;
    printf("<pile2CNV> Bin coverage: binSize=%s | minCoverage=%s%s\n", binSize, minCov, legacyBins ? " | legacy bin labels" : "") > "/dev/stderr";
    printf("<pile2CNV> heteroSNP: minDepth=%s | %s <= normalVAF <= %s\n", minDepth, minVAF, maxVAF) > "/dev/stderr";
    if (snpFile !~ /\.snp/) snpFile = snpFile ".snp";
    printf("<pile2CNV> Writing snpData to %s\n", snpFile) > "/dev/stderr";
    samples = (NF - 3 - hasExonPos) / 2;
    print("<pile2CNV> ", samples, "samples detected") > "/dev/stderr";
    baseHeader = (keepExonPos) ? "Chr\tPos\tExonPos" : "Chr\tPos";
    printf("%s", baseHeader);
    printf("%s", baseHeader) > snpFile;
    for (s = 0; s++ < samples;) {
        printf("\tCov%s", s);
        printf("\tDepth%s\tVAF%s", s, s) > snpFile;
        COL[s] = 2 + (2 * s);
        COL["Read" s] = 3 + (2 * s);
        COVSUM[s] = 0;
    }
    printf("\n");
    printf("\n") > snpFile;
    haveBin = 0;
    next;
}

############# het SNPs: rows where the normal (sample 1) has a non-ref read (SNV, insertion I/i, deletion D/d) #####
$5 ~ /[AaCcTtGgDdIi]/ {
    isGood = 1;
    for (s = 0; s++ < samples;) {
        alt = gsub(/[AaCcTtGgDdIi]/, "", $COL["Read" s]);
        DEPTH[s] = $COL[s];
        if (DEPTH[s] < minDepth) { isGood = 0; break; }
        VAF[s] = (DEPTH[s]) ? alt / DEPTH[s] : -1;
        if ((minVAF > VAF[1]) || (VAF[1] > maxVAF)) { isGood = 0; break; }
    }
    if (isGood) {
        outstring = (keepExonPos) ? $1 "\t" $2 "\t" $NF : $1 "\t" $2;
        for (s = 0; s++ < samples;) outstring = outstring "\t" DEPTH[s] "\t" VAF[s];
        print(outstring) > snpFile;
    }
}

############# coverage bins #############
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
    binChr = $1;
    bin = thisBin;
    # legacy: ExonPos extrapolated to the bin start; fixed: ExonPos of the first covered base (monotonic)
    if (keepExonPos) binExon = legacyBins ? $NF - (pos % binSize) : $NF;
    for (s = 0; s++ < samples;) COVSUM[s] = $COL[s];
    haveBin = 1;
}
{
    pos = $2;
    thisBin = int(pos / binSize);
    if (!haveBin) {
        startBin();
    } else if ($1 != binChr) {
        # chromosome change: the last bin of a chromosome always carries its own start
        flushBin(bin * binSize);
        startBin();
    } else if (thisBin != bin) {
        flushBin(legacyBins ? thisBin * binSize : bin * binSize);
        startBin();
    } else {
        for (s = 0; s++ < samples;) COVSUM[s] += $COL[s];
    }
}
END {
    if (haveBin) flushBin(bin * binSize);
}'
