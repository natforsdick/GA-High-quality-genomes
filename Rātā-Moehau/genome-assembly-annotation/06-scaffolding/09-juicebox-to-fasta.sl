#!/bin/bash

#SBATCH --job-name=JBATout-fasta 
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

# Converting curated Juicebox output back to FASTA

##########
# PARAMS #
INDIR=/PATH/TO/OUTPUTS/output/scaffolding/yahs/
REFDIR=/PATH/TO/OUTPUTS/output/scaffolding/
REF=genome.nextpolish2.fa
YAHSJUICE=/PATH/TO/yahs/juicer

##########
# MODULES
ml purge && ml Python/3.10.5-gimkl-2022a LASTZ/1.04.03-GCC-9.2.0

##########
cd $INDIR
# juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa
$YAHSJUICE post -o genome.nextpolish2-mapped.PT_JBAT-out-postsynteny2 genome.nextpolish2-mapped.PT_JBAT-postsynteny2.review.assembly genome.nextpolish2-mapped.PT_JBAT.liftover.agp ${REFDIR}${REF}
