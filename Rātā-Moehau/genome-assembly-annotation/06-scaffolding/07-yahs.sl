#!/bin/bash

#SBATCH --job-name=rata-yahs 
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G 
#SBATCH --time=00:15:00 
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

##########
# Passing aligned Omni-C data to YAHS to scaffold the assembly

##########
# PARAMS
YAHS='/PATH/TO/yahs/yahs'
IN_DIR='/PATH/TO/OUTPUTS/output/scaffolding/'
IN_BAM='genome.nextpolish2-mapped.PT'
OUT_DIR='/PATH/TO/OUTPUTS/output/scaffolding/yahs/'
REF_DIR='/PATH/TO/OUTPUTS/output/scaffolding/'
REF='genome.nextpolish2.fa'

##########
if [ ! -e ${OUT_DIR} ]; then
	mkdir -p ${OUT_DIR}
else
	echo Found ${OUT_DIR}
fi
cd ${OUT_DIR}

##########
echo Starting YAHS for ${IN_BAM} to scaffold ${REF}
$YAHS ${REF_DIR}${REF} ${IN_DIR}${IN_BAM}.bam -o ${IN_BAM} --no-mem-check
echo Completed YAHS scaffolding
