#!/bin/bash -e

#SBATCH -J STAR-align
#SBATCH --cpus-per-task=16
#SBATCH --mem=14G
#SBATCH -t 00:45:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# aligning RNA-seq data to soft-masked assembly for transcriptome annotation
# based on https://github.com/kherronism/rewarewaannotation

INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/
INDEX=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/alignments/STAR-indexes
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/alignments/STAR-aligned/
PREFIX=metBart-contam-excl.clean
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
REF=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta

# for the guide on readgroup info: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
ml purge; ml STAR/2.7.10b-GCC-11.3.0-alpha

#echo beginning STAR aligment
STAR \
        --genomeDir $INDEX \
        --readFilesIn ${INDIR}sortmerna-out/out/Rata01AG1183009-rrnafilt_R1.fq,${INDIR}sortmerna-out2/out/Rata02AG1183010-rrnafilt_R1.fq \
        ${INDIR}sortmerna-out/out/Rata01AG1183009-rrnafilt_R2.fq,${INDIR}sortmerna-out2/out/Rata02AG1183010-rrnafilt_R2.fq \
        --runThreadN $SLURM_CPUS_PER_TASK \
        --outFileNamePrefix ${OUTDIR}${PREFIX}. \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline ID:Rata01AG1183009 PL:ILLUMINA SM:Rata LB:Rata01AG1183009 , ID:Rata02AG1183010 PL:ILLUMINA SM:Rata LB:Rata02AG1183010 \
        --outSAMtype BAM Unsorted

echo mapped

echo tidying up
if [ -f ${OUTDIR}${PREFIX}.Unmapped.out.mate1 ]; then
        mv ${OUTDIR}${PREFIX}.Unmapped.out.mate1 ${outdir}${prefix}.unmapped_1.fastq
        gzip ${OUTDIR}${PREFIX}.unmapped_1.fastq
fi
if [ -f ${OUTDIR}${PREFIX}.Unmapped.out.mate2 ]; then
        mv ${OUTDIR}${PREFIX}.Unmapped.out.mate2 ${outdir}${prefix}.unmapped_2.fastq
        gzip ${OUTDIR}${PREFIX}.unmapped_2.fastq
fi
echo mapping completed

ml purge
ml SAMtools/1.19-GCC-12.3.0
echo sorting BAM output
samtools sort -o ${OUTDIR}${PREFIX}.Aligned.sorted.bam -T ${OUTDIR}${PREFIX}.sorting --threads 12 ${OUTDIR}${PREFIX}.Aligned.out.bam
echo sorted

# using Picard to collect mapping stats
ml purge
ml picard/2.26.10-Java-11.0.4
echo collecting mapping stats
picard \
        CollectAlignmentSummaryMetrics \
        -I ${OUTDIR}${PREFIX}.Aligned.sorted.bam \
        -O ${OUTDIR}${PREFIX}-metrics.txt \
        -R ${REFDIR}${REF}
echo mapping pipeline complete
