#!/bin/bash -e

#SBATCH -J trimgalore
#SBATCH -c 4
#SBATCH --mem=1G
#SBATCH --array=0-18#
#SBATCH --time=06:00:00 
#SBATCH --output=%x.%j.%a.out
#SBATCH --output=%x.%j.%a.err

# Running trimgalore for adapter trimming and quality filtering
# PARAMS
INDIR=/PATH/TO/rata-pop-gen/data/input/
OUTDIR=/PATH/TO/rata-pop-gen/data/output/02-trimmed/
SAMPLE_LIST=($(<${INDIR}R1-list.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

base=$(basename $SAMPLE _R1_001.fastq.gz) 

ml purge
ml TrimGalore/0.6.10-gimkl-2022a-Python-3.11.3-Perl-5.34.1

cd $INDIR

trim_galore --paired --length 50 \
--three_prime_clip_R1 5 --three_prime_clip_R2 5 --clip_R1 20 --clip_R2 20 --2colour 20 \
--fastqc ${INDIR}${base}_R1_001.fastq.gz ${INDIR}${base}_R2_001.fastq.gz \
-o ${OUTDIR} #--basename ${base}
# trims paired end reads of a sample to a minimum length of 50, does a 3' clip of 5 bp and 5' clip of 20 bp,
# performs clips with a 2-colour compatible quality phred score of 20, and clip of 20 bases

wait
echo "done trimming ${base}"
