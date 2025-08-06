#!/bin/bash 

ml purge
ml BEDTools/2.30.0-GCC-11.3.0

REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/
REF=metBart-contam-excl

cd /PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/

# create masked genome FASTA files
# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi ${REFDIR}${REF}.fasta -bed ${REF}.simple_mask.gff3 \
  -fo ${REF}.simple_mask.soft.fasta

# create complex repeat hard-masked genome
bedtools maskfasta -fi ${REF}.simple_mask.soft.fasta \
  -bed ${REF}.complex_mask.gff3 \
  -fo ${REF}.simple_mask.soft.complex_mask.hard.fasta

# create full soft-masked version - this is what we will use as reference input for gene annotation
bedtools maskfasta -soft -fi ${REF}.simple_mask.soft.fasta \
  -bed ${REF}.complex_mask.gff3 \
  -fo ${REF}.simple_mask.soft.complex_mask.soft.fasta
