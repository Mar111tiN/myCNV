#!/bin/sh

# 1 CORE

### cleans samtools mpileup output 
mawk ' # cat - <(echo "STOP\n") $SNP |   ###### if the SNP markers are necessary, they have to be read in here 
## $5 ~ /[.,]*[ACTG]+[.,]*/ | mawk    # if you only want to output mutant positions (complete the quotes!!)
NR == 1 { 
    minVAF="'${1-0}'";
    ###### QUERY ##############
    # get the letters to look for
    len = split("Aa,Gg,Cc,Tt,Dd,Ii",Letters,",")
    ###### HEADER ################
    for (l=0;l++<len;) {
        LETPAT[l] = "[" Letters[l] "]"
    }
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "Chr","Start","Ref", "Depth", "Alt", "VAF");
    ####### HAVE TO ADJUST ###########
    # loop through the letters
}
{   ######### LINES #############
    # loop through the letters
    for (l = 0;l++< len;) {
        COUNT[l] = gsub(LETPAT[l], "", $5);
    }
  ######### OUTPUT #############
  # loop through the letters
    maxcount = 0;
    base = "";
    # first line extra for pretty
    for (l=0;l++<len;) {
        if (COUNT[l] > maxcount) {
            maxcount = COUNT[l];
            base = substr(Letters[l],1,1);
            vaf =maxcount/$4;
        } 
    }
    if (vaf > minVAF) {
        printf("%s\t%s\t%s\t%s\t",$1,$2,$3,$4);
        printf("%s%s\t%s\n", base,maxcount,vaf);
    } 
}'