#!/bin/sh

minCov=${minCov-0};
keepExonPos=1;





mawk '
#############################################
############# BEGIN #########################
NR == 1 {  ### GET/WRITE HEADER
    ##### ARGS
    binSize='${binSize-10}';
    minCov='${minCov-0}';
    # get the coverage threshold for the output function
    if (minCov) covThresh=minCov * binSize - 1;
    printf("<pile2CNV> Bin coverage: binSize=%s | minCoverage=%s\n", binSize, '$minCov') >> "/dev/stderr";

    # detect XPos
    hasExonPos = ($NF == "ExonPos");
    # keepExonPos only if Exon positions are there
    keepExonPos = (hasExonPos && '${keepExonPos-0}')

    ###### HEADER
    baseHeader = (keepExonPos) ? "Chr\tStart\tExonPos" : "Chr\tStart"
    printf(baseHeader);

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
    print("<pile2CNV> ", samples, "samples detected") > "/dev/stderr"
    # output data col header and init the COVSUM array and create zero-string
    zero="0";
    for (s=0; s++ < samples;) {
        printf("\tCov%s", s);
        # make zero-string
        if (s > 1) {
            zero = zero "|0";
        }
        COVSUM[s] = 0;
    }
    zero=zero "=" zero;

    ### NEWLINE #####
    printf("\n");
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
############# COV DATA ######################
$DEPTHCOL != zero {   ####### COV detection ###############
    pos=$2;
    # split the data to POS and NEG Array
    split($DEPTHCOL, SPLIT, "=");
    split(SPLIT[1], POS, "|");
    split(SPLIT[2], NEG, "|");
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
            }
        }
        lastBin=thisBin;
        
    }
    ####### DEBUG ########
    # print("pos, exonPos, thisBin", pos, exonPos, thisBin)
    ####### DEBUG ########
}
$DEPTHCOL == zero {
    if (!minCount) next;
    print("ZERO");
    printf("%s\t%s", $1, binPos);
    if (keepExonPos) {
        printf("\t%s", exonPos);
        exonPos=$NF-(pos%binSize);
    }
    for (s=0; s++ < samples;) {
        printf("\t0");
    }
    # NEWLINE
    printf("\n");    
}
END {
    outputCov();
}
'