#!/bin/bash

# converting masked outputs to GFF3 
toGFF=/PATH/TO/GFF-CONVERSION-SCRIPT/
cd /PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
REF=metBart-contam-excl

# use Daren's custom script to convert .out to .gff3 for all repeats, simple repeats only, and complex repeats only
# https://github.com/darencard/GenomeAnnotation/blob/master/rmOutToGFF3custom
${toGFF}rmOutToGFF3custom -o ${REF}.full_mask.out > ${REF}.full_mask.gff3
${toGFF}rmOutToGFF3custom -o ${REF}.simple_mask.out > ${REF}.simple_mask.gff3
${toGFF}rmOutToGFF3custom -o ${REF}.complex_mask.out > ${REF}.complex_mask.gff3
