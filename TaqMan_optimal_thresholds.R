##################################################################################################
### Script uses the model underlying eDNAssay to calculate optimal class assignment thresholds 
### for a range of false negative to false positive cost ratios
##################################################################################################
library(caret)

load("TaqMan_trained_model.RData")

thresh <-
  thresholder(trainmodel_taqman, threshold = seq(0, 1, by = 0.01))
colnames(thresh)[5] <- "PPV"
colnames(thresh)[6] <- "NPV"
thresh$PPV[is.nan(thresh$PPV)] <- 1
thresh$FNrate <- 1 - thresh$NPV
thresh$FPrate <- 1 - thresh$PPV

cost <- function(FNrate, FPrate, ratio) {
  print(FNrate * ratio + FPrate * 1)
}

min_cost <- function(FNrate, FPrate, ratio) {
  print(which(cost(thresh$FNrate, thresh$FPrate, ratio) ==
                min(
                  cost(thresh$FNrate, thresh$FPrate, ratio)
                )) / 100)
}

opt_threshold <- c()
for (i in 1:100) {
  opt_threshold[i] <- min_cost(thresh$FNrate, thresh$FPrate, i)
}

save(opt_threshold, file = "TaqMan_optimal_thresholds.RData")
