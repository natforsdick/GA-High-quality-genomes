#!/bin/bash -e

#SBATCH -J filtering
#SBATCH -c 2 #2 cpus
#SBATCH --mem 10G #10 G
#SBATCH -t 10:00:00 # 10 hours
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

# generate various filtered VCFs for comparisons

#########
# PARAMS
INDIR=/PATH/TO/rata-pop-gen/data/output/05-variant-calling/rata-moehau-only/ #directory where files to filter are
vcf_out=${INDIR}filter-trial/
noLD=${vcf_out}noLD/
LD=${vcf_out}LD-filter/
INBCF=rata-moehau-only_VariantCalls_concat_sort.bcf
STATS=/PATH/TO/rata-pop-gen/data/output/05-variant-calling/rata-moehau-only/filter-trial/stats/

#########
ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

#########
cd $INDIR
mkdir -p $vcf_out $noLD $LD $STATS

#for loop to filter file with different values for parameters including
#missingness, depth, and GQ

base=$(basename ${INBCF} _concat_sort.bcf)

# can adjust other filters e.g., maxDP, minQ, as appropriate for the data - best practice is to collect stats from input BCF and assess and adjust from there
#for i in {3..5} #filtering files for 3x, 4x, and 5x depth,  .. means all numbers between those two numbers
for i in {10,15,20}
do
    echo "Filtering SNPs for ${base}...." 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_maf0.01.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.01 \
        --minQ 20 \
        --remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.1site_missing_maf0.01.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.9 \
        --maf 0.01 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_maf0.0.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.0 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.1site_missing_maf0.0.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.9 \
        --maf 0.0 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_maf0.25.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.25 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
done
wait

ml BCFtools/1.10.2-GCC-9.2.0

echo "Filtering for Linkage parameters..."
#for loop to filter previous filtered files for linkage
for bcf in ${noLD}*.bcf
do
    base=$(basename ${bcf} .bcf)
    echo "Running light LD pruning at 0.8 for ${base}...."
    bcftools +prune \
    -l 0.8 \
    -w 1000 \
    -O b \
    -o ${LD}${base}_0.8LD_VariantCalls.bcf \
    ${bcf} 
done
wait

ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

#calculating statistics for no linkage filtered files
for file in ${noLD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --site-depth 
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --depth 
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --missing-site 
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --missing-indv 
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --het
done
wait

#calculating statistics for linkage filtered files
for file in ${LD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --site-depth 
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --depth 
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --missing-site 
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --missing-indv 
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATS}${base} \
        --het
done
wait
