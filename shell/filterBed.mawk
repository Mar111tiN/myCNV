#!/bin/sh

# file with genomic pos in Pos | filterBed <bedfile> [chrom ..which chrom to use] [=0|1..useExonicCoords]
# filters any position-based file to positions included in a bed file
# takes bedfile and as parameter
# works only on stdin in a pipe

bedFile=$1;
echo "STOPBED" > stop.file

cat ${bedFile} stop.file - | mawk '
BEGIN {
    ## INIT ######
    # useExonicCoords 
    useExonicCoords='${3-0}';
    if (useExonicCoords == 1) {
        printf("Printing out exonic coordinates\n") > "/dev/stderr";
    }
    # get filter chrom as arg3 for matching in filterBed file
    # "7" --> chr7
    # default ""  --> chr
    filterChrom="'${2-""}'";
    if (filterChrom !~ "chr") {
        filterChrom = "chr" filterChrom;
    }
    if (filterChrom != "chr"){
        printf("Filtering chromosome %s\n", filterChrom) > "/dev/stderr";
    }
    readBed=1;
    bedCount=0; 
    bedPointer=1;
    bedStep=1;
    chromCount = 1;
}

readBed {
    if ($0 !~ "STOPBED") { # reading bedFile line with right chromosome
        # store the bed regions as blocks in BEDSTART and BEDEND
        if ( $1 !~ filterChrom ) {
            next;
        }
        bedCount++;
        # next chrom in case filterChrom is "chr"
        if ($1 != currentChrom) {
            exonicCoord=0; # the current accumulated exon base count 
            CHROMEND[currentChrom] = bedCount - 1 # fix the last bedPointer for each chrom
            currentChrom = $1;
            CHROMSTART[currentChrom] = bedCount; # for every newly encountered chromosome, create a marker for the bedPointer
            CHROMNAME[bedCount] = currentChrom; # for that position, remember the chrom name
            CHROMCOUNT[chromCount++] = currentChrom;
        }
        
        BEDSTART[bedCount] = $2;
        BEDEND[bedCount] = $3;
        # exonicCoord will be increased by the length of the bed fragment
        # BEDSUM marks the exonic coords of the bed fragments end (defined by bedCount)
        exonicCoord = exonicCoord + $3 - $2;
        # #####
        # print(bedCount,BEDSTART[bedCount], BEDEND[bedCount], exonicCoord, BEDSUM[bedCount]);
        # #####

        BEDSUM[bedCount] = exonicCoord;
    } else { # reached end of bedfile
            # switch to HEADER mode
            readBed = 0;
            writeHeader = 1;
    }
    next;
}

########## HEADER #######################
writeHeader { # check for existence of header and print out
  CHROMEND[currentChrom] = bedCount
  writeHeader = 0;
  readData = 1;
  chromCount = 1;
  currentChrom = CHROMNAME[1]; # reset the current chrom to the first one

  # for (cc in CHROMSTART) {
  #   print(cc, CHROMSTART[cc]);
  # }
    if ($0 ~ "Chr\t") { 
        printf($0);
        if (useExonicCoords == 1) {
            printf("\tExonPos");
        }
            printf("\n");
            next;
        } else { # move on to readData on same line
            print("<filterBedAll> No header detected") > "/dev/stderr";
        }
}

############# DATA #########################
readData {  # switching to data
    # get data
    chrom=$1;
    pos=$2;
    
    # check if chrom is contained in bedfile
    if (chrom in CHROMSTART == 0) next;
    # move to the right chromosome in the read file
    # should best be already done in the samtools file
    while (chrom != currentChrom) {
        currentChrom = CHROMCOUNT[++chromeCount];
        bedPointer = CHROMSTART[currentChrom];
    }
  
    # cycle through the bedRegions
    while (pos >= BEDSTART[bedPointer] ) { # if pos downstream of current bedRegion, drop line silently
        #print(bedPointer, BEDSTART[bedPointer], BEDEND[bedPointer], pos);
        if (pos > BEDEND[bedPointer]) { # step to next bedRegion
            bedPointer++;
            # if all bedRegions have been read, stop output
            if (bedPointer > bedCount) exit 0;
            # if bedpointer is in next chrom, skip
            if (bedPointer > CHROMEND[currentChrom]) next;
        } else { # write to file
            # DEBUG############
            # print(bedPointer, BEDSUM[bedPointer], BEDSTART[bedPointer], BEDEND[bedPointer]);
            # DEBUG ########## 
            # print the base data
            printf($0);
            if (useExonicCoords==1) {
                exonicCoord=BEDSUM[bedPointer] - BEDEND[bedPointer] + pos;
                printf("\t%s", exonicCoord);
            }
            printf("\n"); 
            break; 
        }
    }
}'
# cleanup
# cd rm stop.file;