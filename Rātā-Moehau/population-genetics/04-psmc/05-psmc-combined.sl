#!/bin/bash -e

#SBATCH -J psmc-combine
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH -t 00:15:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# Define output directory
OUTDIR=/PATH/TO/rata-pop-gen/data/psmc/

# Chage to output directory
cd $OUTDIR

# load the psmc module
module purge; module load psmc/0.6.5-gimkl-2018b

# Define the base name (without file extension)
base_name="EXT049-08_S8_map_rata_tahae_ref_diploid"

# Combine the bootstrap results into one file
echo "Combining results for $base_name"
cat ${base_name}-v3.psmc ${base_name}_round-*.psmc > ${base_name}_combined-v3.psmc

# Plot using the combined file with the updated settings
echo "Plotting combined results"
#psmc_plot.pl -u 2.06e-9 -g 30 -R ${base_name}-v3-boot -png ${base_name}_combined-v3.psmc

#brassicaceae -g 20
#psmc_plot.pl -u 1.304e-7 -g 20 -R ${base_name}-v3-boot-g20-brass ${base_name}_combined-v3.psmc

#walnut -g 20
psmc_plot.pl -u 4.12e-8 -g 20 -R ${base_name}-v3-boot-g20-wal ${base_name}_combined-v3.psmc

#arabidposis -g 20
#psmc_plot.pl -u 7.0e-9 -g 20 -R ${base_name}-v3-boot-g20-arab ${base_name}_combined-v3.psmc

#prunus -g 20
#psmc_plot.pl -u 9.5e-9 -g 20 -R ${base_name}-v3-boot-g20-prun ${base_name}_combined-v3.psmc  

#echo "Plotting single results"
# Plot using the original single file with the updated settings
#psmc_plot.pl -u 2.06e-9 -g 30 -R ${base_name}-single ${base_name}-v3.psmc

