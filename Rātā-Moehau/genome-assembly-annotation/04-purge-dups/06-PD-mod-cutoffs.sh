#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA, modified by Nat Forsdick

# step 06: modify cutoffs

##########
# PARAMS
PURGE_DUPS=/PATH/TO/purge_dups/bin/
OUTDIR=/PATH/TO/OUTPUTS/output/04-purge-dups/asm-shasta/
PRE=asm-shasta # PREFIX

R1=01-
R2=03- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
CUTOFFS="-l 35 -m 94 -u 228" # -u 240" # as determined by considering coverage plot produced in previous step
##########

cd ${OUTDIR}
echo $CUTOFFS
${PURGE_DUPS}calcuts ${CUTOFFS} ${R1}${PRE}-PB.stat > ${R2}${PRE}-cutoffs

# Following this, you need to run steps 04-07 with $ROUND modified for new cutoffs.
