#!/bin/sh

########### GENOMIC COVERAGE ############
# takes normal bam files downstream of samtools view 
# processes the cigar string for coverage computation
# PARAMETER:
#     $1 is minimal coverage to be output
# USAGE:
# samtools bam_file [region] | bamCoverage 0
# OUTPUT:
#   Chr Pos Coverage
# !!! only outputs positions where the coverage is changing
# this makes it sparse

mawk '

BEGIN {
  minCov="'${1-0}'"; # get $1 with 0 as default
  cigPat = "^[0-9]+[NMDIS]";  
  # init pointer
  pos = 0;
  cov = 0;
  # print Header
  printf("%s\t%s\t%s\n","Chr","Pos","Coverage");
}

# look at all hits
$6 !~ /\*/ {
  # get the data
  chrom = $3;
  start = $4; # start position is pushed up by cigar info
  cigar = $6;

  readPos = start # readPos is current read position during cigar parse

  ############ START - END ############################
  # get the start end end coords for the coverage from the cigar
  while (match(cigar, cigPat)) {
    # length of cigar block
    l = substr(cigar,RSTART,RLENGTH-1)
    # type of cigar block
    t = substr(cigar,RSTART+RLENGTH-1,1)
    # if match block, push coords into coverage array
    if (t == "M") {
      # with every read, total coverage increases by one between start and end
      COVup[readPos] ++;
      readPos += l; # readPos is moved up by the length of the Mcigar
      COVdown[readPos] ++; # coverage goes down at new start pos
    # jump down over intron gaps
    } else if (t ~ /[ND]/) {
      readPos += l; # N and D result in no coverage
    }
    # reduce cigar string
    cigar = substr(cigar, RSTART+RLENGTH);
  }

  ########### OUTPUT ################
  ########### SPARSE ################
  if (pos < start) {
    for (p=pos; p++< start - 1;) { # walk over last coords
      change = COVup[p] - COVdown[p];
      cov += change;
      if (cov >= minCov && change != 0) {  # this makes it sparse output
        printf("%s\t%s\t%s\n", chrom, p,cov);
      } 
      delete COVup[p];
      delete COVdown[p];
    }
    pos=start-1; # go up to start -1 as other reads at that position can increase
    # .. coverage
    stop=readPos;
  }
}
END { # spill out all the remaining data arrays
  for (p=pos; p++ < stop;){
    change = COVup[p] - COVdown[p];
    cov += change;
    if (cov >= minCov && change != 0) {  # this makes it sparse output
      printf("%s\t%s\t%s\n", chrom, p,cov);
    } 
  }

}
'