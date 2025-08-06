#!/bin/bash -e

#SBATCH -J plot-psmc
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G 
#SBATCH -t 0:30:00 
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

OUTDIR=/PATH/TO/rata-pop-gen/data/psmc
cd $OUTDIR

module purge; module load psmc/0.6.5-gimkl-2018b
echo plotting
# different versions of the .psmc files originate from setting different parameter space in the previous script
# comment out as desired
#psmc_plot.pl -u 2.06e-9 -g 30 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v1-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v1.psmc
#psmc_plot.pl -u 2.06e-9 -g 30 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v2-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v2.psmc
#psmc_plot.pl -u 2.06e-9 -g 30 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 2.06e-9 -g 30 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot EXT049-08_S8_map_rata_tahae_ref_diploid_combined-v3.psmc
#psmc_plot.pl -u 2.06e-9 -g 30 -pY50000 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot-exp EXT049-08_S8_map_rata_tahae_ref_diploid_combined-v3.psmc
#psmc_plot.pl -u 2.06e-9 -g 10 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot-g10 EXT049-08_S8_map_rata_tahae_ref_diploid_combined-v3.psmc

# scaling the mutation rate per generation, using the rate from the walnut paper
#psmc_plot.pl -u 2.06e-9 -g 1 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g1-u2.06e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 4.12e-9 -g 2 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g2-u4.12e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 1.03e-8 -g 5 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g5-u1.03e-8-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 2.06e-8 -g 10 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g10-u2.06e-8-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 4.12e-8 -g 20 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u4.12e-8-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 6.18e-8 -g 30 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g30-u6.18e-8-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 2.06e-7 -g 100 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g100-u2.06e-7-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc

# trying high and low mutationrates from the angiosperm vs gymnosperm paper, this is from brassicaceae, highest rate
#psmc_plot.pl -u 6.52e-9 -g 1 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g1-u6.52e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc

# trying to match what Choi et all used for their metrosideros paper which includes demographic modelling.
#psmc_plot.pl -u 7.0e-9 -g 20 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u7.0e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 7.0e-9 -g 25 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g25-u7.0e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 9.5e-9 -g 25 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g25-u9.5e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
#psmc_plot.pl -u 9.5e-9 -g 20 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u9.5e-9-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc

# getting some other generation times for brassicaceae, which has a higher mutation rate compared to walnut
psmc_plot.pl -u 1.304e-7 -g 20 -R EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u1.3e-7-plot EXT049-08_S8_map_rata_tahae_ref_diploid-v3.psmc
