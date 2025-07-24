#!/bin/bash -e 

#SBATCH -J omnic_qc 
#SBATCH -c 2 
#SBATCH --mem=2G 
#SBATCH --time=00:05:00 
#SBATCH --output %x.%j.out #  
#SBATCH --error %x.%j.err # 

###############
# Run Omni-C QC using the bam file produced from mapping Omni-C reads to assembly.
# This script takes two arguments, $1 : the path to the input directory, / 
# and $2 : the prefix of the input bam file. 
# e.g. to execute: 
# sbatch 02-omnic-qc.sl /PATH/TO/OUTPUTS/rata-omnic-qc/ Rata-QC_clean1_R-sorted 

############### 
module purge
module load Miniconda3/4.10.3
ml SAMtools/1.13-GCC-9.2.0
 
###############
# ENVIRONMENT 
source /PATH/TO/Miniconda3/4.10.3/etc/profile.d/conda.sh
conda activate hic_qc
hic_qc=/PATH/TO/hic_qc/hic_qc.py
IN_DIR=$1
IN_BAM=$2
CPU=20

#####################
cd ${IN_DIR}
# Check that sorted BAM is present
if [ ! -e ${IN_DIR}${IN_BAM}-sorted.bam ]; then
	echo "Sorting ${HIC}"
	samtools sort -n -@ $CPU ${IN_DIR}${IN_BAM}.bam -o ${IN_DIR}${IN_BAM}-sorted.bam
else
	echo "BAM found, running QC"
fi

###############
cd $IN_DIR
echo "running hic_qc for ${IN_BAM}"
python ${hic_qc} -b ${IN_BAM}-sorted.bam -o ${IN_BAM}-sorted.hicqc 
echo "finished running hic_qc for ${IN_BAM}"
