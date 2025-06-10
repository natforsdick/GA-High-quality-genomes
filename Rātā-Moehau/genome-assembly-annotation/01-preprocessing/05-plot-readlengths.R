# Rscript to plot the outputs of 04-get-readlengths.sh
sessionInfo()
library('dplyr')
library(kableExtra)

Run1 <- read.table("/PATH/TO/INPUT/output/rata-MinION/Rata_1/sup-fastq/combined-sup-fastqs/rata-PB5-pass_read_length.txt")
Run2 <- read.table("/PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/rata-PB5-pass_read_length.txt")

# Rename columns
dfList <- list(Run1,Run2)
dflist <- lapply(dfList, function(x) {
  names(x)[ grep("V2", names(x))] <- "Length"
  names(x)[ grep("V1", names(x))] <- "Count"
  x} )

library(kableExtra)
# Get summary for each data set in list
summary1 <- lapply(dflist,summary)
summary1 %>% 
  kbl() %>%
  kable_styling()

Run1 <- rename(Run1, Length = V2, Count = V1)
Run2 <- rename(Run2, Length = V2, Count = V1)

# Add dataset name as column
Run1 <- data.frame(append(Run1, c(Date='2022-06-11')))
Run2 <- data.frame(append(Run2, c(Date='2023-03-20')))

plot(Run1[1:5000,],type="l")

library(ggplot2)

myplots<-lapply(dflist,function(x) 
  p <- ggplot(x,aes(Length, Count)) + geom_line(col="blue")+
    xlab("Length (bp)") +
    ylab("Count")
)

myplots

RataLength <- rbind(Run1, Run2)

g <- ggplot(RataLength, aes(Length, Count))
g + geom_line(aes(col=factor(Date)), alpha=0.8) + 
  labs(x="Length (bp)",
       col="Batch")

# looking more closely at the short end of the spectrum:
g + geom_line(aes(col=factor(Date)), alpha=0.8) + 
  labs(x="Length (bp)",
       col="Batch") +
  xlim(0,25000)
