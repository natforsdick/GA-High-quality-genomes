#!/bin/bash -e

#SBATCH -J bwa
#SBATCH --time=02:30:00 
#SBATCH --mem=28G 
#SBATCH --cpus-per-task=32
#SBATCH --array=1-16
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

# An array to iterate the mapping process over all samples

###########
# PARAMS
# need to make new fq_list for trimgalore outputs
fq_list=/PATH/TO/rata-pop-gen/data/output/03-merged/EXT049P2/merged-samplist.txt

reffile=metBart-contam-excl
refdir=/PATH/TO/rata-pop-gen/data/genome/
ref=$refdir$reffile

INDIR=/PATH/TO/rata-pop-gen/data/output/03-merged/EXT049P2/
SAMDIR=/PATH/TO/rata-pop-gen/data/output/04-mapped/EXT049P2/sam/
BAMDIR=/PATH/TO/rata-pop-gen/data/output/04-mapped/EXT049P2/bam/
fq1=_val_1.fq.gz #Read 1 suffix
fq2=_val_2.fq.gz #Read 2 suffix

###########
# MODULES
module purge
module load  BWA/0.7.17-GCC-9.2.0 SAMtools/1.10-GCC-9.2.0
###########

# MAPPING
cd $INDIR

QUERY1=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${fq_list} ) #`cat ${fq_list} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $1}'`

echo "Task ${SLURM_ARRAY_TASK_ID}; Processing $QUERY1"

# capture readgroup info
base=$(basename $QUERY1 _val_1.fq.gz)
infoline=$(zcat ${QUERY1} | head -n 1)
instrument=`echo ${infoline} | cut -d ':' -f1`
instrumentrun=`echo $infoline | cut -d ':' -f2`
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`
platform="Illumina"

# incorporate sample information into the alignment
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${base}_library1"
rgsm="SM:${base}"

echo "Aligning reads for $base"
bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 32 $ref $QUERY1 ${datadir}${base}${fq2} > ${SAMDIR}${base}.sam

echo "Converting sam file to bam file for $base"
samtools view -@ 32 -T $ref.fa -b ${SAMDIR}${base}.sam > ${BAMDIR}${base}.bam

echo "Sorting and indexing file"
samtools sort -@ 32 -o ${BAMDIR}${base}.aligned.sorted.bam ${BAMDIR}${base}.bam

echo "Removing intermediate files"
rm ${SAMDIR}${base}.sam
rm ${BAMDIR}${base}.bam
echo "Completed $QUERY1"
