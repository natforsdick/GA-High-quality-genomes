#!/bin/bash -e

#SBATCH -J braker-rna
#SBATCH --cpus-per-task=18
#SBATCH --mem=32G
#SBATCH -t 6:30:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# Running BRAKER3 for rata Moehau annotation using aligned RNAseq data as evidence
# does require an Augustus config file to be prepared

INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/alignments/STAR-aligned/
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/braker/rna/
PREFIX=metBart-contam-excl
REF=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
BRAKER=/PATH/TO/BRAKER/3.0.8-gimkl-2022a-Perl-5.34.1/scripts/
BAMIN=${INDIR}${PREFIX}.clean.Aligned.sorted.bam # input BAM alignment of cleaned RNAseq data

ml purge; ml BRAKER/3.0.8-gimkl-2022a-Perl-5.34.1 compleasm/0.2.5-gimkl-2022a

export PYTHONPATH=/PATH/TO/compleasm/0.2.5-gimkl-2022a:$PYTHONPATH

cd $OUTDIR
braker.pl \
        --genome ${REFDIR}${REF} \
        --species metBart \
        --workingdir $OUTDIR \
        --bam $BAMIN \
        --softmasking_off \
        --threads 48 \
        --busco_lineage=eudicots_odb10 \
        --AUGUSTUS_SCRIPTS_PATH=/PATH/TO/AUGUSTUS/3.5.0-gimkl-2022a/scripts \
        --AUGUSTUS_BIN_PATH=/PATH/TO/AUGUSTUS/3.5.0-gimkl-2022a/bin \
        --AUGUSTUS_CONFIG_PATH=${OUTDIR}config
