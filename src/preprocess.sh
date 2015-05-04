#!/bin/sh
# $1 threads
# $2 tmpdir
# $3 fasta
# $4 bams

# run mpileup
echo "mpileup"
cut -f1 $3.fai | parallel -j $1 "samtools mpileup -g -u -t DP -f $3 -r {} $4 | bcftools annotate -O u --remove 'INFO/INDEL,INFO/IDV,INFO/IMF,INFO/I16,INFO/QS' - > $2/{}.bcf"
# run concat
echo "concat"
cut -f1 $3.fai | xargs -I "{}" bcftools concat -O b "$2/{}.bcf"
# delete tmpfiles
cut -f1 $3.fai | xargs -I "{}" rm "$2/{}.bcf"
