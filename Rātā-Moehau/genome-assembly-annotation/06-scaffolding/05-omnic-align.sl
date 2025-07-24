#!/bin/bash -e

#SBATCH --job-name=rata-omnic-map
#SBATCH --cpus-per-task=18
#SBATCH --mem=12G 
#SBATCH --time=15:00:00 
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

##########
# MODULES
ml purge
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0 pairtools/1.0.2-gimkl-2022a-Python-3.10.5

##########
# PARAMS
PREFIX=genome.nextpolish2
INDIR='/PATH/TO/OUTPUTS/output/scaffolding/' 
OMNICR1=/PATH/TO/OUTPUTS/output/scaffolding/rata-omnic-clean-R1.fastq.gz
OMNICR2=/PATH/TO/OUTPUTS/output/scaffolding/rata-omnic-clean-R2.fastq.gz
REF=genome.nextpolish2
TMPDIR="/PATH/TO/tmp-omnic-${SLURM_JOB_ID}"
CPU=18
##########

cd $INDIR
mkdir $TMPDIR
export TMPDIR=$TMPDIR

# alignment -T0 = all alignments to generate stats, do downstream quality filtering
echo aligning
bwa mem -5SP -T0 -t $CPU $REF.fa $OMNICR1 $OMNICR2 -o ${PREFIX}-aligned.sam

##########
# find ligation junctions
echo finding ligation junctions
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $CPU \
--nproc-out $CPU --chroms-path ${REF}.genome ${PREFIX}-aligned.sam > ${PREFIX}-parsed.pairsam

##########
# sort pairsam
echo sorting pairsam
pairtools sort --nproc $CPU --tmpdir=$TMPDIR ${PREFIX}-parsed.pairsam > ${PREFIX}-sorted.pairsam

##########
# remove duplicates
echo removing duplicates
pairtools dedup --nproc-in $CPU --nproc-out $CPU --mark-dups --output-stats ${PREFIX}-stats.txt \
--output ${PREFIX}-dedup.pairsam ${PREFIX}-sorted.pairsam

##########
# split .bam, .pairs
echo splitting bam
pairtools split --nproc-in $CPU --nproc-out $CPU --output-pairs ${PREFIX}-mapped.pairs \
--output-sam ${PREFIX}-unsorted.bam ${PREFIX}-dedup.pairsam

##########
if [ -f ${PREFIX}-unsorted.bam ]
then
echo "pipeline completed"
else
echo "pipeline not complete"
fi
