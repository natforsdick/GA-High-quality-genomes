#!/bin/bash -e

#SBATCH --job-name=fastp 
#SBATCH --cpus-per-task=12
#SBATCH --mem=4G 
#SBATCH --time=02:10:00 
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

# use Fastp to trim and quality filter the raw Omni-C data from the full sequencing run

##########
# PARAMS
OUTDIR=/PATH/TO/OUTPUTS/output/scaffolding/
ASSEMBLY=/PATH/TO/OUTPUTS/output/polish/NextPolish/genome.nextpolish2.fa
APREFIX=Rata-OmniC-full # output prefix
HIC_DIR=/PATH/TO/INPUT/rata-omnic/
HIC_RAW1=${HIC_DIR}FILENAME_S1_L001_R
HIC_RAW2=${HIC_DIR}FILENAME_S1_L002_R
READ1=1_001.fastq.gz
READ2=2_001.fastq.gz
##########

ml purge && module load fastp/0.23.2-GCC-11.3.0

##########
# Clean HiC Reads with fastp
echo processing ${HIC_RAW1}$READ1
fastp \
-i ${HIC_RAW1}${READ1} \
-o ${HIC_DIR}${APREFIX}1_clean1_R1.fastq.gz \
-I ${HIC_RAW1}${READ2} \
-O ${HIC_DIR}${APREFIX}1_clean1_R2.fastq.gz \
--trim_front1 15 \
--trim_front2 15 \
--qualified_quality_phred 20 \
--length_required 50 \
--thread 12
#--trim_tail1 15 \ # only for MiSeq QC data
#--trim_tail2 15 \
#--cut_tail 

# uncomment the following if Omni-C data generated via NovaSeq where you get two sets of reads back
#echo processing ${HIC_RAW2}$READ1
#fastp \
#-i ${HIC_RAW2}${READ1} \
#-o ${HIC_DIR}${APREFIX}2_clean1_R1.fastq.gz \
#-I ${HIC_RAW2}${READ2} \
#-O ${HIC_DIR}${APREFIX}2_clean1_R2.fastq.gz \
#--trim_front1 15 \
#--trim_front2 15 \
#--qualified_quality_phred 20 \
#--length_required 50 \
#--thread 12
