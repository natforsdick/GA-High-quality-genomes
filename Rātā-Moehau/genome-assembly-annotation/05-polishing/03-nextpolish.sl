#!/bin/bash -e

#SBATCH --job-name nextpolish
#SBATCH --cpus-per-task=20
#SBATCH --mem=36G
#SBATCH --time=03:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

##########
# PARAMS #
INDIR=/PATH/TO/PROCESSED-ILLUMINA/output/illumina/02-trimmed/
R1=${INDIR}EXT041-06_S6-all-R1.fq.gz
R2=${INDIR}EXT041-06_S6-all-R2.fq.gz
OUTDIR=/PATH/TO/OUTPUTS/output/polish/NextPolish/
REFDIR=/PATH/TO/MEDAKA-OUT/output/polish/medaka/
REF=${REFDIR}consensus.fasta # round 1
round=2
threads=20
NEXTPOLISH=/PATH/TO/NextPolish/lib/nextpolish1.py
#REF=genome.nextpolish1.fa # round 2
##########

#mkdir -p $OUTDIR # round 2 comment out
cd $OUTDIR

ml purge && ml Miniconda3 SAMtools/1.13-GCC-9.2.0 BWA/0.7.17-gimkl-2017a
source /PATH/TO/Miniconda3/4.8.3/etc/profile.d/conda.sh
conda activate NextPolish # see whether it can use paralleltask in the nextpolish workflow

for ((i=1; i<=${round};i++)); do
#step 1:
   echo index the genome file and do alignment $i;
   bwa index ${REF};
   bwa mem -t ${threads} ${REF} ${R1} ${R2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3 - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam;
   echo index bam and genome files $i;
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${REF};
   echo polish genome file $i;
   python ${NEXTPOLISH} -g ${REF} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp${i}.fa;
   echo round $i complete;
#step2:
   echo index genome file and do alignment $i;
   bwa index $REF;
   bwa mem -t ${threads} ${REF} ${R1} ${R2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam;
   echo index bam and genome files $i;
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${REF};
   echo polish genome file $i;
   python ${NEXTPOLISH} -g ${REF} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish${i}.fa;
   echo round $i complete;
done;
#Finally polished genome file: genome.nextpolish2.fa
