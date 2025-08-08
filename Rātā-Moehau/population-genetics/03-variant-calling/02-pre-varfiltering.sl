#!/bin/bash -e

#SBATCH -J pre-filt
#SBATCH -c 16
#SBATCH --mem=4G
#SBATCH --time=2:30:00 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Prepping vcfs for filtering
bcfdir=/PATH/TO/rata-pop-gen/data/output/05-variant-calling/rata-moehau-only/ #bcf file directory
samplist=/PATH/TO/rata-pop-gen/data/output/04-mapped/rata-moehau-only-bams/bam-samplist.txt

ml purge
ml BCFtools/1.15.1-GCC-11.3.0

for file in ${bcfdir}*.vcf
do
	base=$(basename $file .vcf)
	echo Processing ${file}

	#echo compressing ${file}
        bgzip -c ${bcfdir}${base}.vcf > ${bcfdir}${base}.bcf.gz
	echo indexing ${file}
	tabix ${bcfdir}${base}.bcf.gz
	# put bcf files names into a list for concatenation
	ls ${bcfdir}${base}.bcf.gz >> ${bcfdir}bcflist.txt
done

# concatenate the chunked bcf files into a whole population bcf
# -a allow overlaps, -D remove exact duplicates (outputs a single record for any duplicates
echo concatenating
bcftools concat --file-list ${bcfdir}bcflist.txt -a -D -O b -o ${bcfdir}rata-moehau-only_VariantCalls_concat.bcf --threads 24

# now let's sort the resultant bcf
ml purge; ml BCFtools/1.19-GCC-11.3.0
TMPDIR=${bcfdir}temp
mkdir -p $TMPDIR

bcftools sort --temp-dir $TMPDIR -O b -o ${bcfdir}rata-moehau-only_VariantCalls_concat_sort.bcf ${bcfdir}rata-moehau-only_VariantCalls_concat.bcf
