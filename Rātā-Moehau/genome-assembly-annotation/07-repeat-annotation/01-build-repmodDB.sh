#!/bin/bash

# build the repeat modeler database for the assembly

REFDIR=/PATH/TO/OUTPUTS/output/scaffolding/
REF=metBart-contam-excl.fasta 
species=metBart
OUTDIR=/PATH/TO/OUTPUTS/07-annotation/repeats/

mkdir -p $OUTDIR
cd $OUTDIR

ml purge && ml RepeatModeler/2.0.3-Miniconda3

BuildDatabase -name $species $REFDIR$REF
