#!/bin/bash -e

#SBATCH --job-name=rata-omnic-sort
#SBATCH --cpus-per-task=10 
#SBATCH --mem=30G 
#SBATCH --time=00:10:00
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err

# sort the aligned Omni-C data ready for scaffolding

##########
# MODULES
ml purge && ml SAMtools/1.15.1-GCC-11.3.0

##########
# PARAMS
PREFIX=genome.nextpolish2
INDIR='/PATH/TO/OUTPUTS/output/scaffolding/'
OMNICR1=/PATH/TO/OUTPUTS/output/scaffolding/rata-omnic-clean-R1.fastq.gz
OMNICR2=/PATH/TO/OUTPUTS/output/scaffolding/rata-omnic-clean-R2.fastq.gz
REF=genome.nextpolish2
TMPDIR="/PATH/TO/tmp-omnic-${SLURM_JOB_ID}"
CPU=20

##########
cd $INDIR
mkdir $TMPDIR
export TMPDIR

# sort bam
echo sorting bam
samtools sort -@${CPU} -T ${TMPDIR}tempfile.bam -o ${PREFIX}-mapped.PT.bam ${PREFIX}-unsorted.bam

##########
# index bam
echo indexing final bam
samtools index ${PREFIX}-mapped.PT.bam
echo omnic processing complete
