#!/bin/bash -e

# Prepare the assembly file for alignment
# Following a first round of assembly scaffolding with YaHS and subsequent manual curation, 
# you may wish to align the Omni-C data again for a second round of scaffolding, in which
# case you will need to modify the params.

ml purge  && ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0

INDIR=/PATH/TO/OUTPUTS/output/polish/NextPolish/
REF=genome.nextpolish2 # assembly prefix
OUTDIR=/PATH/TO/OUTPUTS/output/scaffolding/

mkdir -p $OUTDIR
cd $OUTDIR

# symlink the assembly to the outdir for convenience
ln -s ${INDIR}${REF}.fa ${OUTDIR}${REF}.fa

# make samtools index
samtools faidx $REF.fa

# make .genome file
cut -f1,2 $REF.fa.fai > $REF.genome

# make BWA index
bwa index $REF.fa
