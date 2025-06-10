#!/bin/bash

# get readlengths for the combined raw ONT sequence data

PASSDIR=/PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/

cd $PASSDIR
zcat rata-PB5-pass.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > rata-PB5-pass_read_length.txt
