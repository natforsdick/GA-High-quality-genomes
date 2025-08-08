#!/bin/bash -e

#SBATCH -J mpileup
#SBATCH -c 16
#SBATCH --mem=6G
#SBATCH --time=1-12:00:00 #Walltime (HH:MM:SS) 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Run BCFTOOLS mpileup tool for variant calling
ref=/PATH/TO/rata-pop-gen/data/genome/metBart-contam-excl.fa
bamdir=/PATH/TO/rata-pop-gen/data/output/04-mapped/rata-moehau-only-bams/
bcfdir=/PATH/TO/rata-pop-gen/data/output/05-variant-calling/rata-moehau-only/ # output bcf file directory
samplist=/PATH/TO/rata-pop-gen/data/output/04-mapped/rata-moehau-only-bams/bam-samplist.txt

ml purge; ml BCFtools/1.15.1-GCC-11.3.0 SAMtools/1.15.1-GCC-11.3.0

cd $bamdir 

if [ ! -e ${samplist} ]; then
	ls -d ${bamdir}*.bam > samplist.txt
fi

if [ ! -e ${bcfdir}chunks/ ]; then
	mkdir -p ${bcfdir}chunks
fi

# chunk bam files into 6 pieces using custom perl script 
# @Lanilen/SubSampler_SNPcaller/split_bamfiles_tasks.pl
# chunked files will help mpileup run faster
echo chunking
perl /PATH/TO/rata-pop-gen/scripts/03-variant-calling/split_bamfiles_tasks.pl -b ${samplist} -g $ref -n 6 -o ${bcfdir}chunks | parallel -j 6 {}

#echo chunking complete
#run bcftools mpileup in parallel on chunks of bam files with BCFtools
echo running mpileup
for (( i=1; i<=6; i++ )); do
        bcftools mpileup -E -O b -f $ref -a AD,ADF,DP,ADR,SP -o ${bcfdir}rata_${i}_raw.bcf ${bcfdir}chunks/${i}/* &
done
wait
echo mpileup complete

#SNP calling on bcf files with bcftools call
for file in ${bcfdir}*.bcf
	do
	base=$(basename $file .bcf)
	bcftools call $file --threads 24 -mv -O v -o ${bcfdir}${base}_VariantCalls.vcf &
done
wait
echo variant calling is complete
