#!/bin/bash -e

#SBATCH --job-name=make-paf
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:45:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=32

# Purge_dups pipeline
# Created by Sarah Bailey, UoA, modified by Nat Forsdick

# step 01: align HiFi sequencing data to the assembly and generate a paf file
# Takes one parameter - PRI or ALT (primary or alternate haplotype)

#########
# MODULES
module purge
module load minimap2/2.24-GCC-11.3.0

#########
# PARAMS
INDIR=/PATH/TO/INPUT/output/03-assembly/asm-shasta/
OUTDIR=/PATH/TO/OUTPUT/output/04-purge-dups/asm-shasta/
DATA=/PATH/TO/PROCESSED/ONT/output/rata-MinION/
PRE=asm-shasta # assembly prefix
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-

cd $OUTDIR

echo "Indexing"
date
if [ ! -e ${INDIR}${PRE}.mmi ]; 
then
minimap2 -d ${INDIR}${PRE}.mmi ${INDIR}${PRE}.fasta
else
echo "index found"
fi

echo "Mapping"
date
minimap2 -x map-ont -t 24 ${INDIR}${PRE}.mmi \
${DATA}rata-1-trimmed.fastq.gz ${DATA}rata-2-trimmed.fastq.gz | gzip -c - > ${R1}${PRE}-mapped.paf.gz
echo "done"
date
