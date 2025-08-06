#!/bin/bash -e

#SBATCH -J repmod
#SBATCH --cpus-per-task=28
#SBATCH --mem=12G
#SBATCH -t 1-02:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# running RepeatModeler using the database generated from the assembly to identify de novo repeats

##########
# PARAMS
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/
species=metBart

##########
# MODULES
ml purge && ml RepeatModeler/2.0.3-Miniconda3

cd $OUTDIR
# -pa = parellel search - each pa uses 4 cpu
RepeatModeler -database $species -pa 24 -LTRStruct
