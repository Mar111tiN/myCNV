#!/bin/sh

### cleans output from samtools mpileup (for one sample) and only returns Chr, Start, and the read data
# Chr   Start   Depth Read1
# this can be done with pon2cols tool or deliberate cutting
# removes from all reads the traces of base location and indel lengths

mawk '
{
    read = $5;
    # remove position traces from all read fields
    gsub(/\^[^\t]|\$/,"",read);
    # 
    while (match(read,/[+-][0-9]+/)) {
        pre = substr(read,1,RSTART-2);
        indel = substr(read,RSTART,1); # + or -
        base = substr(read,RSTART-1,1); # base before the deletion
        l = substr(read,RSTART+1,RLENGTH-1);
        post = substr(read,RSTART+RLENGTH+l);
        if (indel == "-") {
            if (match(base,/[ACGT]/)) {
                base = "D";
            } else {
                base = "d";
            }
        } else {
            if (match(base,/[ACGT]/)) {
                base = "I";
            } else {
                base = "i";
            }            
        }

        read = pre base post;
    }     
# print all fields
    for (i=0;i++<4;) {
        printf("%s\t",$i);
    };
    printf("%s\n",$5);
}' | mawk '$5 ~ /[.,]*[ACTG]+[.,]*/'   # only output actual SNPs
