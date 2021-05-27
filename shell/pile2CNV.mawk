#!/bin/sh
# v1.0
# computes coverages and SNP output from pileup for multiple samples

# INPUT
# takes clean pileup from samtools mpileup with ref (...,,,AT.Ii..)
# and base qualities removed
# can contain ExonPos (auto-detected)
# 1 <= samples < N
# 	Chr	Start	Ref	Cov1	Read1	Cov2	Read2	...	CovN	ReadN	[ ExonPos ]	
# first sample should be the normal sample for SNP VAF

# USAGE: 
# samtools mpileup -f $HG38 -q 20 -Q 25 -l $BED -r chr? | cleanpileup | filterBed $BED -x -c chr? | pile2CNV
# [     -m | --min-coverage]        <INT=0>                     minimum (binned) coverage for output    ]
# [     -x | --keep-exon-pos        <Flag=False>                if exonPos should be used               ]
# [     -b | --bin-size             <INT=10>                    size of the bin for summed coverage     ]
# [     -o | --output-snp-file      <path to snp file>          snp_file for the heteroSNP output       ]
# [     -v|--min-vaf|--min-VAF      <INT=0>                     minimum normal VAF for heteroSNP output ]
# [     -V|--max-vaf|--max-VAF      <FLOAT=0.75>                maximum normal VAF for heteroSNP output ]
# [     -d|--min-depth              <INT=10>                    minimum depth (of any sample) heteroSNP ]                

# OUTPUT
#   SNP
#       Chr Start [ExonPos] Ref Depth1 VAF1 Depth2	VAF2...
#   COV
#       Chr Start [ExonPos] Cov1 Cov2


####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--keep-exon-pos)
        keepExonPos=1
        shift
        ;;
        # snp_file output
        -o|--output-snp-file)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            snpFile=$2
            shift 2
        else
            echo "<pile2CNV> Error: output file for heteroSNP is missing\n[-o|--output-snp-file]" >&2
            exit 1
        fi
        ;;
        # minimum coverage
        -c|--min-coverage)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minCov=$2
            shift 2
        else
            echo "<pile2CNV> Error: Argument for $1 is missing" >&2
            exit 1
        fi
        ;;
        # bin Size
        -b|--bin-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            binSize=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide size for step size [-s|--step-size (default=10)]" >&2
            exit 1
        fi
        ;;
        # minVAF
        -v|--min-vaf|--min-VAF)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minVAF=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide minimum VAF for heteroSNP output [-v|--min-vaf (default=0)]" >&2
            exit 1
        fi
        ;;
        # maxVAF
        -V|--max-vaf|--max-VAF)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            maxVAF=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide maximum VAF for heteroSNP output [-v|--min-vaf (default=0)]" >&2
            exit 1
        fi
        ;;
        # minDepth
        -d|--min-depth)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minDepth=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide minimum VAF for heteroSNP output [-d|--min-depth (default=10)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<pile2CNV> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done








################ MAWK #########################
###############################################
mawk '
############# BEGIN #########################
NR == 1 {  ### GET/WRITE HEADER
    ##### ARGS
    # detect ExonPos
    hasExonPos = ($NF == "ExonPos");
    # keepExonPos only if Exon positions are there
    keepExonPos = (hasExonPos && '${keepExonPos-0}')

    # COV
    minCov='${minCov-0}';
    binSize='${binSize-10}';
    # get the coverage threshold for the output function
    if (minCov) covThresh=minCov * binSize - 1;
    printf("<pile2CNV> Bin coverage: binSize=%s | minCoverage=%s\n", binSize, '$minCov') >> "/dev/stderr";

    # SNP
    minVAF='${minVAF-0.1}';
    maxVAF='${maxVAF-0.9}';
    minDepth='${minDepth-10}'
    printf("<pile2CNV> heteroSNP: minDepth=%s | %s <= normalVAF <= %s\n", minDepth, minVAF, maxVAF) >> "/dev/stderr";
    snpFile="'${snpFile-"test.snp"}'";
    if (snpFile !~ ".snp") snpFile = snpFile ".snp"
    # setup SNP file and write its Header
    printf("<pile2CNV> Writing snpData to %s\n", snpFile) >> "/dev/stderr";

    ###### INIT #########
    # detect samples
    samples = (NF-3-hasExonPos)/2;
    print("<pile2CNV> ", samples, "samples detected") > "/dev/stderr"


    ###### HEADER
    # print headers to stream and to SNP output file
    baseHeader = (keepExonPos) ? "Chr\tStart\tExonPos" : "Chr\tStart"

    printf(baseHeader);
    printf(baseHeader) > snpFile;

    # output data col header
    for (s=0; s++ < samples;) {
        printf("\tCov%s", s);
        printf("\tDepth%s\tVAF%s", s, s) >> snpFile;
        # COL stores the col number for each sample
        COL[s] = 2 + (2 * s);
        COL["Read" s] = 3 + (2 * s);
        COVSUM[s] = 0;
    }
    ### NEWLINE #####
    printf("\n");
    printf("\n") >> snpFile;
    next;
}

############# VAFDATA #######################
$5 ~ /[AaCcTtGgDdI]/ {  # only look at rows where the normal has a mutation

    ##### DEBUG #########
    # print($0) >> snpFile;
    # print($NormalCol) >> snpFile;
    # print("minDepth", minDepth, "minVAF", minVAF, "maxVAF", maxVAF) >> snpFile;
    # print("samples", samples) >> snpFile;
    # for (s=0;s++<samples;) print(s, $COL[s]) >> snpFile;;
    ##### DEBUG #########

    # isGood as marker for doing output
    isGood=1
    # get the Data
    for (s=0;s++<samples;) {
        alt=gsub(/[AaCcTtGgDdI]/, "", $COL["Read" s]);
        DEPTH[s]=$COL[s];
        # if depth too low (for any sample) break loop and mark lowDepth
        if (DEPTH[s] < minDepth) {
            isGood=0;
            break;
        }
        VAF[s]=(DEPTH[s]) ? alt/DEPTH[s] : -1;
        # filter normal VAFs for minVAF and maxVAF
        #  varcan needs normal_bam tumor_bam input so normal has s==1
        if ((minVAF > VAF[1]) || (VAF[1] > maxVAF)) {
            isGood=0;
            break;
        }
    }
    ## OUTPUT if good
    if (isGood) {
        # build up the output string
        outstring=(keepExonPos) ? $1 "\t" $2  "\t" $NF : $1 "\t" $2;

        for (s=0;s++<samples;) {
            outstring = outstring "\t" DEPTH[s] "\t" VAF[s];
            # printf("\t%s\t%.3f", DEPTH[s], VAF[s]) >> snpFile;
        }
        print(outstring) >> snpFile;
    }
}

#############################################
############# COV DATA ######################
function outputCov(  s) {
    # get positions for thisBin
    binPos=thisBin * binSize;
    # check minimum coverage
    if (minCov) {
        noOutput=1;
        for (s=0; s++ < samples;) {
            if (COVSUM[s] > covThresh) {
                noOutput=0;
                break;
            }
        }

        if (noOutput) { # do the update without print
        ### DEBUG
        # print("noOutput");
        ### DEBUG
            if (keepExonPos) exonPos=$NF-(pos%binSize);
            for (s=0; s++ < samples;) {
                # reset COVSUM
                COVSUM[s] = $COL[s];
            }
            return;
        }
    }
    printf("%s\t%s", $1, binPos);

    if (keepExonPos) {
        printf("\t%s", exonPos);
        exonPos=$NF-(pos%binSize);
    }
    for (s=0; s++ < samples;) {
        printf("\t%s", COVSUM[s] / binSize);
        # reset COVSUM
        COVSUM[s] = $COL[s];
    }
    # NEWLINE
    printf("\n");
}

{   ####### COV detection ###############
    pos=$2;
    # get the current bin
    thisBin=int((pos) / binSize);

    if (thisBin == lastBin) { # same bin --> accumulate
        for (s=0; s++ < samples;) {
            COVSUM[s] += $COL[s];
        }
    } else { # jump to next bin
        ########### OUTPUT LAST BIN AND UPDATE ##########
        if (lastBin > 0) {
            outputCov()
        } else { # get the positions for first bin
            if (keepExonPos) {
                binPos=thisBin * binSize;
                exonPos=$NF-(pos%binSize);
            }
        }
        lastBin=thisBin;
    }
    ####### DEBUG ########
    # print("pos, exonPos, thisBin", pos, exonPos, thisBin)
    ####### DEBUG ########
}
END {
    outputCov();
}
'
