#!/bin/sh
# v1.0.0
# computes coverages and SNP output from PON matrix

# INPUT
# takes PONmatrix with optional ExonPos added by filterBED (auto-detected)
# 	Chr	Start	Ref	A	G	C	T	I	D	Depth	[ ExonPos ]	

# USAGE: 
# samtools gunzip <  chr1.pon.gz [ | filterBed $BED -x -c chr? ] | PON2CNV
# [     -c | --min-coverage]        <INT=0>                     minimum (binned) coverage for output    ]
# [     -x | --keep-exon-pos        <Flag=False>                if exonPos should be used               ]
# [     -b | --bin-size             <INT=10>                    size of the bin for summed coverage     ]
# [     -o | --output-snp-file      <path to snp file>          snp_file for the heteroSNP output       ]
# [     -v|--min-vaf|--min-VAF      <INT=0>                     minimum total VAF for heteroSNP output  ]
# [     -d|--min-depth              <INT=10>                    minimum total depth for heteroSNP output]                

# OUTPUT
#   SNP
#       Chr Pos [ExonPos] Depth VAF               Depth and VAF are for aggregated over all samples
#   COV
#       Chr Pos [ExonPos] Cov1 Cov2..CovN


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
            echo "<PON2CNV> Error: output file for heteroSNP is missing\n[-o|--output-snp-file]" >&2
            exit 1
        fi
        ;;
        # minimum coverage
        -c|--min-coverage)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minCov=$2
            shift 2
        else
            echo "<PON2CNV> Error: Argument for $1 is missing" >&2
            exit 1
        fi
        ;;
        # bin Size
        -b|--bin-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            binSize=$2
            shift 2
        else
            echo "<PON2CNV> Error: Provide size for step size [-s|--step-size (default=10)]" >&2
            exit 1
        fi
        ;;
        # minVAF
        -v|--min-vaf|--min-VAF)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minVAF=$2
            shift 2
        else
            echo "<PON2CNV> Error: Provide minimum VAF for heteroSNP output [-v|--min-vaf (default=0.25)]" >&2
            exit 1
        fi
        ;;
        # minDepth
        -d|--min-depth)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minDepth=$2
            shift 2
        else
            echo "<PON2CNV> Error: Provide minimum VAF for heteroSNP output [-d|--min-depth (default=50)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<PON2CNV> Error: Unsupported flag $1" >&2
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
    printf("<PON2CNV> Bin coverage: binSize=%s | minCoverage=%s\n", binSize, '$minCov') >> "/dev/stderr";

    # SNP
    VAFPRECISION=4;
    PREC=10^VAFPRECISION;
    minVAF='${minVAF-0.25}';
    minDepth='${minDepth-50}'
    printf("<PON2CNV> heteroSNP: minDepth=%s | totalVAF <= %s\n", minDepth, minVAF) >> "/dev/stderr";
    snpFile="'${snpFile-"test.snp"}'";
    if (snpFile !~ ".snp") snpFile = snpFile ".snp"
    # setup SNP file and write its Header
    printf("<PON2CNV> Writing snpData to %s\n", snpFile) >> "/dev/stderr";


    ###### HEADER
    baseHeader = (keepExonPos) ? "Chr\tPos\tExonPos" : "Chr\tPos"
    printf(baseHeader);
    printf(baseHeader) > snpFile;
    
    printf("\tVAF\tDepth\n") >> snpFile;
    # detect DEPTH column in fields
    for (f=0; f++<NF;) {
        if ($f == "Depth") DEPTHCOL=f;
    }
    # reset lastBin;
    lastBin=0;

    next;
}

NR == 2 { # detect the sample count in the first data line and output the data HEADER
    split($DEPTHCOL, SPLIT, "=");
    samples = gsub(/\|/, "", SPLIT[1]) + 1;
    print("<PON2CNV> ", samples, "samples detected") > "/dev/stderr"
    # output data col header and create zero-string
    zero="0";
    
    for (s=0; s++ < samples;) {
        printf("\tCov%s", s);
        # make zero-string
        if (s > 1) {
            zero = zero "|0";
        }
        # init POS, NEG arrays for use as global vars in function
        POS[s]=0;
        NEG[s]=0;
    }
    zero=zero "=" zero;
    ### NEWLINE #####
    printf("\n");

    ### VAF #######
    totalDEPTH=0;
}

function aggString(col,  SPLIT, sum) {
    # aggregate all values for one data col
    sum=0;
    split(col, SPLIT, "=");
    split(SPLIT[1], POS, "|");
    split(SPLIT[2], NEG, "|");
    for (s=0; s++<samples;) {
        sum+=POS[s] + NEG[s];
    }
    return sum;
}


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


#############################################
############# DATA ######################
$DEPTHCOL != zero { # skip absolute zero lines
    pos=$2;

    ##### VAF ##################
    # get the total ALT values
    totalAlt=0;
    for (col=3; ++col<DEPTHCOL;) {
        # skip zero strings
        if ($col == zero) continue;
        # get the split data arrays
        totalAlt += aggString($col);
    }
    # get the total Depths and fill POS and NEG with Depth values
    # ..for use in coverage 
    totalDepth=aggString($DEPTHCOL)
    totalVAF=int((totalAlt / totalDepth) * PREC) / PREC;
    if ((totalVAF > minVAF) && (totalDepth >= minDepth)) {
        # build up the output string
        outstring=$1 "\t" pos;
        if (keepExonPos) outstring = outstring "\t" $NF;
        outstring = outstring "\t" totalVAF "\t" totalDepth;
        print(outstring) >> snpFile;
    }

    ##### COV #################
    # get the current bin
    thisBin=int((pos) / binSize);

    if (thisBin == lastBin) { # same bin --> accumulate
        for (s=0; s++ < samples;) {
            COVSUM[s] += POS[s] + NEG[s];
        }
    } else { # jump to next bin
        ########### OUTPUT LAST BIN AND UPDATE ##########
        if (lastBin > 0) {
            outputCov()
        } else { # first line
            # get the positions for first bin
            binPos=thisBin * binSize;
            if (keepExonPos) { 
                exonPos=$NF-(pos%binSize);
            }
            # store first line
            for (s=0; s++ < samples;) {
                COVSUM[s] = POS[s] + NEG[s];
                totalDepth += COVSUM[s];
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