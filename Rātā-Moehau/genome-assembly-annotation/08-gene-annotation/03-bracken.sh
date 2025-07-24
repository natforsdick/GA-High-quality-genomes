#!/bin/bash

# Run Braken to interpret Kraken reports

##########
# PARAMS
INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/02-kraken2/
samples="Rata01AG1183009 Rata02AG1183010"
readlength=150
k2db=/PATH/TO/KRAKENDB/databases/k2_standard_20240605

##########
ml purge; ml Kraken2/2.1.3-GCC-11.3.0 Bracken/2.7-GCC-11.3.0

cd $INDIR

# -l classification level, default 'S' = species.
# -t threshold, default 10. For species with <=10 reads will not receive any additional reads from higher taxonomy levels when distributing reads for abundance estimation
for samp in $samples
do
    bracken -d $k2db -i ${samp}.kreport -o ${samp}.bracken -r ${readlength} -l S -t 10
done
