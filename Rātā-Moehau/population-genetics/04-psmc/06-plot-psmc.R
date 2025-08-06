## How to plot PSMC results
# Following this tutorial: https://psmc-tutorial-birdlab.readthedocs.io/en/latest/03_psmc_visualization.html

#load necessary libraries
library("ggpubr")

#select a font for your plots
op <- par(family = "serif")

#lets make a blank plot
# Set scipen to a high value to prevent scientific notation
options(scipen = 999)
plot(1,1,
     ylim=c(0,100),xlim=c(2,150000),
     log="x",type="n",
     main="",
     ylab = expression("Effective Population Size" ~ "   x 10"^4),
     xlab = expression("Years ago"),
     axes = FALSE)  # Turn off default axes

# Add custom x-axis with whole numbers
axis(1, at = c(10, 100, 1000, 10000, 100000), labels = c("10 K", "100 K", "1 M", "10 M", "100 M"))
# Add default y-axis
axis(2)
# Optionally, add a box around the plot
box()

#blank plot with smaller effective population size axis (for showing 4 mutation rates):
options(scipen = 999)
plot(1,1,
     ylim=c(0,25),xlim=c(2,150000),
     log="x",type="n",
     main="",
     ylab = expression("Effective Population Size" ~ "   x 10"^4),
     xlab = expression("Years ago"),
     axes = FALSE)  # Turn off default axes

# Add custom x-axis with whole numbers
axis(1, at = c(10, 100, 1000, 10000, 100000), labels = c("10 K", "100 K", "1 M", "10 M", "100 M"))
# Add default y-axis
axis(2)
# Optionally, add a box around the plot
box()

#blank plot with smaller effective population size axis and shorter y axis (for showing just walnut at g= 20 mutation rates):
options(scipen = 999)
plot(1,1,
     ylim=c(0,5),xlim=c(2,15000),
     log="x",type="n",
     main="",
     ylab = expression("Effective Population Size" ~ "   x 10"^4),
     xlab = expression("Years ago"),
     axes = FALSE)  # Turn off default axes

# Add custom x-axis with whole numbers
axis(1, at = c(10, 100, 1000, 10000), labels = c("10 K", "100 K", "1 M", "10 M"))
# Add default y-axis
axis(2)
# Optionally, add a box around the plot
box()

#blank plot with larger effective population size axis and shorter y axis (for showing just walnut at range of mutation rates):
options(scipen = 999)
plot(1,1,
     ylim=c(0,100),xlim=c(2,15000),
     log="x",type="n",
     main="",
     ylab = expression("Effective Population Size" ~ "   x 10"^4),
     xlab = expression("Years ago"),
     axes = FALSE)  # Turn off default axes

# Add custom x-axis with whole numbers
axis(1, at = c(10, 100, 1000, 10000), labels = c("10 K", "100 K", "1 M", "10 M"))
# Add default y-axis
axis(2)
# Optionally, add a box around the plot
box()


#lets add lines at different time point of our interest.
abline(v=18, col="black")
abline(v=1040, col="black")
#abline(v=c(6,22,120), col="black")
#text(119.85,4.18514,"LIG",cex=1.0)
text(49,4.18514,"Ōtira Glaciation, LGM",cex=1.0)
text(2090,4.58514,"tombolo forms",cex=1.0)
#text(5.928641,4.18514,"MDH",cex=1.0)

# extract bootstrap values for each species separately.this is for generation time of 30 years (incorrectly scaled)
for(i in 0:96) {
  path<-paste0("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot.",i,".txt")
  boot.iter<-read.table(path)
  ya<-boot.iter[-1,1]
  ya.pl=ya[which(ya>4999)]
  ne<-boot.iter[-1,2]
  ne.pl=ne[which(ya>4999)]
  lines((ya/1000),ne,xlog=T,type="l",col=alpha("darkblue",0.1),lwd=0.8,ylim=c(0,50))

# extract bootstrap values for each species separately.this is for generation time of 10 (incorrectly scaled)
for(i in 0:96) {
  path<-paste0("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot-g10.",i,".txt")
  boot.iter<-read.table(path)
  ya<-boot.iter[-1,1]
  ya.pl=ya[which(ya>4999)]
  ne<-boot.iter[-1,2]
  ne.pl=ne[which(ya>4999)]
  lines((ya/1000),ne,xlog=T,type="l",col=alpha("darkgreen",0.1),lwd=0.8,ylim=c(0,50))
}
  

# extract bootstrap values for walnut u for 20 year generation time, correctly scaled i.e. u = 4.12 x 10-8
for(i in 0:96) {
    path<-paste0("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot-g20-wal.",i,".txt")
    boot.iter<-read.table(path)
    ya<-boot.iter[-1,1]
    ya.pl=ya[which(ya>4999)]
    ne<-boot.iter[-1,2]
    ne.pl=ne[which(ya>4999)]
    lines((ya/1000),ne,xlog=T,type="l",col=alpha("grey",0.1),lwd=0.8,ylim=c(0,50))
  }  
  

#adding the legend to the plot - feel free to customize.
legend(150, 85, legend=c( "g=10", "g=30"),
       col=c( "green", "cornflowerblue"), lty=1, cex=1.0, box.lty=5)

# lets extract values of the main PSMC results for each generation time separately.
g30 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-plot.0.txt")
g30.ya=g30[-1,1]
g30.ne=g30[-1,2]

g10 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-boot-g10.0.txt")
g10.ya=g10[-1,1]
g10.ne=g10[-1,2]

#generation = 1, mutation rate = 2.06e-9 Walnut U
g1u2.06 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g1-u2.06e-9-plot.0.txt") 
g1u2.06.ya=g1u2.06[-1,1]
g1u2.06.ne=g1u2.06[-1,2]

#generation = 2, mutation rate = 4.12e-9 Walnut u adapted for longer generation time
g2u4.12 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g2-u4.12e-9-plot.0.txt") 
g2u4.12.ya=g2u4.12[-1,1]
g2u4.12.ne=g2u4.12[-1,2]

#generation time = 5, mutation rate = 1.03e-8
g5u1.03 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g5-u1.03e-8-plot.0.txt") 
g5u1.03.ya=g5u1.03[-1,1]
g5u1.03.ne=g5u1.03[-1,2]

#generation time = 10, mutation rate = 2.06e-8
g10u2.06 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g10-u2.06e-8-plot.0.txt") 
g10u2.06.ya=g10u2.06[-1,1]
g10u2.06.ne=g10u2.06[-1,2]

#generation time = 20, mutation rate = 4.12e-8 Walnut u adapted for longer generation time
g20u4.12 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u4.12e-8-plot.0.txt") 
g20u4.12.ya=g20u4.12[-1,1]
g20u4.12.ne=g20u4.12[-1,2]

#generation time = 30, mutation rate = 6.18e-8
g30u6.18 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g30-u6.18e-8-plot.0.txt") 
g30u6.18.ya=g30u6.18[-1,1]
g30u6.18.ne=g30u6.18[-1,2]

#generation time = 100, mutation rate = 2.06e-7
g100u2.06 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g100-u2.06e-7-plot.0.txt") 
g100u2.06.ya=g100u2.06[-1,1]
g100u2.06.ne=g100u2.06[-1,2]

#generation time = 1, mutation rate = 6.52e-9 (brassicaceae)
g1u6.52 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g1-u6.52e-9-plot.0.txt") 
g1u6.52.ya=g1u6.52[-1,1]
g1u6.52.ne=g1u6.52[-1,2]

#genration time = 20, mutation rate = 1.304×10−7
g20u1.3 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u1.3e-7-plot.0.txt")
g20u1.3.ya=g20u1.3[-1,1]
g20u1.3.ne=g20u1.3[-1,2]

#from choi et al gen time = 20, mutation rate = 7 x 10 -9 (arabidopsis)
g20u7 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u7.0e-9-plot.0.txt")
g20u7.ya=g20u7[-1,1]
g20u7.ne=g20u7[-1,2]

#from choi et al gen time = 25, mutation rate = 7 x 10 -9 (arabidopsis)
g25u7 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g25-u7.0e-9-plot.0.txt")
g25u7.ya=g25u7[-1,1]
g25u7.ne=g25u7[-1,2]

#from choi et al gen time = 20, mutation rate = 9.5 x 10 -9 (prunus)
g20u9.5 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g20-u9.5e-9-plot.0.txt")
g20u9.5.ya=g20u9.5[-1,1]
g20u9.5.ne=g20u9.5[-1,2]

#from choi et al gen time = 25, mutation rate = 9.5 x 10 -9 (prunus)
g25u9.5 = read.table("EXT049-08_S8_map_rata_tahae_ref_diploid-v3-g25-u9.5e-9-plot.0.txt")
g25u9.5.ya=g25u9.5[-1,1]
g25u9.5.ne=g25u9.5[-1,2]

#lets plot those extracted values.
#Walnut not scaled for generation time
lines(x=(g30.ya/1000),y=g30.ne,type="l",col="lightblue",lwd=2.5)

lines(x=(g10.ya/1000),y=g10.ne,type="l",col="green",lwd=2.5)

#walnut corrected
lines(x=(g1u2.06.ya/1000),y=g1u2.06.ne,type="l",col="darkblue",lwd=2.5)
lines(x=(g2u4.12.ya/1000),y=g2u4.12.ne,type="l",col="darkgreen",lwd=2.5)
lines(x=(g5u1.03.ya/1000),y=g5u1.03.ne,type="l",col="lightblue",lwd=2.5)
lines(x=(g10u2.06.ya/1000),y=g10u2.06.ne,type="l",col="lightgreen",lwd=2.5)
lines(x=(g20u4.12.ya/1000),y=g20u4.12.ne,type="l",col="black",lwd=2.5)
lines(x=(g30u6.18.ya/1000),y=g30u6.18.ne,type="l",col="blue",lwd=2.5)
lines(x=(g100u2.06.ya/1000),y=g100u2.06.ne,type="l",col="grey",lwd=2.5)

#brassicaceae
lines(x=(g1u6.52.ya/1000),y=g1u6.52.ne,type="l",col="red",lwd=2.5)
lines(x=(g20u1.3.ya/1000),y=g20u1.3.ne,type="l",col="darkgrey",lwd=2.5)

#Arabidopsis
lines(x=(g20u7.ya/1000),y=g20u7.ne,type="l",col="lightgreen",lwd=2.5)
lines(x=(g25u7.ya/1000),y=g25u7.ne,type="l",col="green",lwd=2.5)

#Prunus
lines(x=(g20u9.5.ya/1000),y=g20u9.5.ne,type="l",col="blue",lwd=2.5)
lines(x=(g25u9.5.ya/1000),y=g25u9.5.ne,type="l",col="darkblue",lwd=2.5)


#adding another legend to the plot
legend(2, 95, legend=c( "g=10, u=2.06e-9 (WRONG u)", "g=30, u=2.06e-9 (WRONG u)", "g=1, u=2.06e-9 (Walnut u)", "g=2, u=4.12e-9", "g=5, u=1.03e-8", "g=10, u=2.06e-8", "g=20, u=4.12e-8", "g=30, u=6.18e-8", "g=100, u=2.06e-7", "g=1, u=6.52e-9 (brassicaceae u - high)"),
       col=c("green", "lightblue", "red", "purple", "yellow", "orange","brown","darkgreen", "grey", "black"), lty=1, cex=1.0, box.lty=5)

#legend for smaller plot. with shorter Y axis and fewer lines in it
legend (2,20,legend=c("G=20 u=1.304×10−7 (Brassicaceae)", "G=20 u=4.12x10-8 (walnut)","G=20 u=9.5x10-9 (Prunus)","G=20 u=7x10-9 (Arabidopsis)"),
        col=c("darkgrey", "lightblue", "blue", "lightgreen"), lty=1, cex=1, box.lty=5)

#legend for walnut plot
legend(2,95,legend=c("g=1, u=2.06e-9", "g=2, u=4.12e-9", "g=5, u=1.03e-8", "g=10, u=2.06e-8", "g=20, u=4.12e-8", "g=30, u=6.18e-8", "g=100, u=2.06e-7"),
        col=c("darkblue", "darkgreen", "lightblue", "lightgreen", "black", "blue", "grey"), lty=1, cex=1, box.lty=5)
