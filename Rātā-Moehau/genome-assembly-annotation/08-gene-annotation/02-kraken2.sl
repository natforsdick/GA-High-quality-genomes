#!/bin/bash -e

#SBATCH -J kraken2
#SBATCH --cpus-per-task=6
#SBATCH --mem=36G
#SBATCH -t 00:10:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# For Kraken2 classification of input RNAseq data to identify and exclude contamination

##########
# PARAMS
INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/01-trim/
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/02-kraken2/

mkdir -p $OUTDIR
cd $OUTDIR

ml purge; ml Kraken2/2.1.3-GCC-11.3.0

echo running kraken for one set of paired-end sequences
kraken2 --db $KRAKEN2_DEFAULT_DB --threads 16 --use-names --report Rata01AG1183009.kreport --output Rata01AG1183009-kraken2-out.txt \
    --classified-out Rata01AG1183009-classified-seqs#.fq --unclassified-out Rata01AG1183009-clean-seqs#.fq \
    --paired ${INDIR}Rata01AG1183009_1_val_1.fq.gz ${INDIR}Rata01AG1183009_2_val_2.fq.gz
    
echo running kraken for second set of paired-end sequences
kraken2 --db $KRAKEN2_DEFAULT_DB --threads 16 --use-names --report Rata02AG1183010.kreport --output Rata02AG1183010-kraken2-out.txt \
    --classified-out Rata02AG1183010-classified-seqs#.fq --unclassified-out Rata02AG1183010-clean-seqs#.fq \
    --paired ${INDIR}Rata02AG1183010_1_val_1.fq.gz ${INDIR}Rata02AG1183010_2_val_2.fq.gz
