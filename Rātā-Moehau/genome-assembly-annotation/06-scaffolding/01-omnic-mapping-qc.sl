#!/bin/bash 

#SBATCH --job-name=qc-mapping 
#SBATCH --cpus-per-task=24
#SBATCH --mem=4G
#SBATCH --time=00:05:00 
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

###############
# Mapping of Omni-C QC data generated from MiSeq, following the Arima Genomics mapping pipeline to map one paired-end dataset
# https://github.com/ArimaGenomics/mapping_pipeline/

###############
# MODULES
###############
module purge
module load BWA/0.7.17-GCC-9.2.0 picard/2.21.8-Java-11.0.4 SAMtools/1.13-GCC-9.2.0 samblaster/0.1.26-GCC-9.2.0

###############
# PARAMS
###############
HIC='Rata-QC_clean1_R' #'basename_of_fastq_files'
LABEL='rata-omnicqc-mapped' #'overall_exp_name'
BWA='bwa' #'/path/to/bwa/bwa-0.7.12/bwa'
SAMTOOLS='samtools' #'/path/to/samtools/samtools-1.3.1/samtools'
IN_DIR='/PATH/TO/INPUT/rata-omnic-qc/' # '/path/to/gzipped/fastq/files'
REF_DIR='/PATH/TO/OUTPUTS/output/polish/NextPolish/'
REF='/PATH/TO/OUTPUTS/output/polish/NextPolish/genome.nextpolish2.fa' #'/path/to/reference_sequences/reference_sequences.fa'
FAIDX='$REF.fai'
PREFIX='genome.nextpolish2' #'bwa_index_name'
RAW_DIR='/PATH/TO/OUTPUTS/rata-omnic-qc/' #'/path/to/write/out/bams'
STATS='/PATH/TO/scripts/Hi-C_scripts/get_stats.pl' #'/path/to/get_stats.pl'
PICARD='/PATH/TO/picard/2.21.8-Java-11.0.4/picard.jar'
TMP_DIR='/PATH/TO/tmp/' #'/path/to/write/out/temporary/files'
MAPQ_FILTER=10
CPU=20
############

echo "Checking output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
# Indexing takes around 40 min for a 1.2 Gb genome
cd $REF_DIR
if [ -f ${REF}.amb ]; then
    echo "${REF} index found"
    else
	bwa index -a bwtsw -p $PREFIX $REF
	echo "Finished indexing $REF"
fi

echo "Starting QC alignment"
echo "bwa mem -t $CPU -5SP $REF ${IN_DIR}${HIC}1.fastq.gz ${IN_DIR}${HIC}2.fastq.gz | samblaster |\
 samtools view -@ $CPU -buSh -F 2316 - > ${RAW_DIR}${HIC}.bam"

bwa mem -t $CPU -5SP $REF ${IN_DIR}${HIC}1.fastq.gz ${IN_DIR}${HIC}2.fastq.gz | samblaster |\
 samtools view -@ $CPU -buSh -F 2316 - > ${RAW_DIR}${HIC}.bam

cd ${RAW_DIR}
echo "Sorting ${HIC}"
samtools sort -@ $CPU ${RAW_DIR}${HIC}.bam -o ${RAW_DIR}${HIC}-sorted.bam
echo "Indexing ${HIC}"
samtools index -@ $CPU ${RAW_DIR}${HIC}-sorted.bam

echo "Finished QC alignment"
