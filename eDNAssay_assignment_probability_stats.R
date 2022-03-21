##################################################################################################
### Script processes specificity data to output summary stats on multiple sequences per taxon
##################################################################################################
library(data.table)
setwd()

### Assignment probability stats
aps <-
  read.csv(file.choose(), header = TRUE) # Specificity file output by eDNAssay
aps <- data.table(aps)

aps[, ':='(
  Min_amp = min(Amp),
  Max_amp = max(Amp),
  Mean_amp = mean(Amp),
  SD_amp = sd(Amp),
  N = length(Amp)
), by = Taxon]
aps[is.na(aps)] <- 0
aps <- aps[order(Taxon),]
aps <- unique(aps, by = "Taxon")
aps <- aps[,-c(3:4)]

write.csv(aps, "Assay_stats.csv", row.names = FALSE)
print("Finished!")
