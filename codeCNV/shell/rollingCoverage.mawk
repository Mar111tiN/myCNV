#!/bin/sh

# ################ GENOMIC ROLLING COVERAGE ##############
# takes the output from chromCoverage and outputs the rolling average coverage per 100 (default) bases


mawk '
NR == 1 {
  width=int("'${1:-100}'" / 2) * 2; # parameter for the rolling scope (default = 100)
  half=width / 2;
  # get the coords for the fields
  # for flexible data structures
  n = split("Chr,Pos,Coverage",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  } 
  # print Header
  printf("%s\t%s\t%s\n","Chr","Pos","Coverage");
}

# special for the first data row
NR == 2 { 
  ##### INIT ###########
  # READ LINE
  lastPos=$COORD["Pos"];
  lastCov=$COORD["Coverage"];
  # posCenter is the middle between sumL and sumR
  posCenter=int(lastPos / half) * half;
}

NR > 2 {
  # READ LINE
  chrom=$COORD["Chr"];
  pos=$COORD["Pos"];
  cov=$COORD["Coverage"];

  # first check, whether a posCenter border has been crossed
  # sumL is left of posCenter [posCenter-half..posCenter)
  # sumR is right of posCenter including posCenter [posCenter..posCenter + half)


  if (pos >= posCenter + half) { # jumped over a border 
      # fill up sumB
      sumR += lastCov * (half - lastPos + posCenter);
      covMean = (sumL + sumR) / width;
      printf("%s\t%s\t%s\n",chrom,posCenter, covMean);

    if (pos >= posCenter + width) { # stepped over two borders
      # spill the right side coverage
      sumL = lastCov * half;
      covMean = (sumR + sumL) / width;
      printf("%s\t%s\t%s\n",chrom,posCenter + half,covMean);
      posCenter=int(pos / half) * half;
  
    } else { # stepped over one border
      # step up the posCenter and fill sumR
      posCenter = posCenter + half;
      sumL = sumR;
    }
    sumR = lastCov * (pos - posCenter + 1);
  } else { # did not jump a border
    sumR += lastCov * (pos - lastPos);
  }
  lastCov = cov;
  lastPos = pos;
  # printf("Pos: %s\tCov: %s\tsumL: %s\tsumR: %s\n",pos, cov, sumL, sumR);
 }
 END {
   print(posCenter,sumL, sumR)
   printf("%s\t%s\t%s\n",chrom,posCenter,(sumL + sumR) / width);
   printf("%s\t%s\t%s\n",chrom,posCenter + half, sumR / width);
 }'