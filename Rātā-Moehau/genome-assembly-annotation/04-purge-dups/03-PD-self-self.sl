#!/bin/bash -e

#SBATCH --job-name=self-self
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=16

# Purge_dups pipeline
# Created by Sarah Bailey, UoA, modified by Nat Forsdick

# step 03: do a self-self alignment
# Takes 1 input parameter: either 'PRI' or 'ALT'

#########
# MODULES
module purge
module load minimap2/2.24-GCC-11.3.0 
#########

#########
# PARAMS
OUTDIR=/PATH/TO/output/04-purge-dups/asm-shasta
PRE=asm-shasta # PREFIX
R1=01- # Designate cutoffs round - either default (01) or modified (02)
R2=02-
#########

cd $OUTDIR
# -x asm5: intra-specific asm-to-asm alignment
if [ "$1" == "PRI" ]; then 
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}${PRI}.split ${R1}${PRE}${PRI}.split | gzip -c - > ${R1}${PRE}${PRI}.split.self.paf.gz
elif [ "$1" == "ALT" ]; then
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}${ALT}.split ${R1}${PRE}${ALT}.split | gzip -c - > ${R1}${PRE}${ALT}.split.self.paf.gz
else
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}.split ${R1}${PRE}.split | gzip -c - > ${R1}${PRE}.split.self.paf.gz
fi 
