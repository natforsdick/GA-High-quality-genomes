#!/bin/bash -e

# index the masked reference assembly file

##########
# PARAMS
INDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
FASTA=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/alignments/STAR-out/

mkdir -p $OUTDIR
cd /PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/alignments/

# calculate scaling
NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${FASTA}.fai`
echo $NUM_BASES

ml purge; ml STAR/2.7.10b-GCC-11.3.0-alpha

echo STAR indexing
# RAM limited to 12 GB
STAR --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles ${INDIR}${FASTA} --runThreadN 4 --genomeSAindexNbases ${NUM_BASES} --limitGenomeGenerateRAM 12000000000
echo indexing complete
