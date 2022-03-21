##################################################################################################
### Script trains a random forest classifier to predict qPCR cross-amplification
##################################################################################################
library(tidyverse)
library(caret)
library(randomForest)
library(beepr)
library(plotROC)

### Select input CSV file containing training data
traindata <-
  read.csv(file.choose(), header = TRUE) # "Training_data"

##################################################################################################
### Define variables of interest and create training dataframe
Assay <- as.factor(traindata$Assay)
Species <- as.factor(traindata$Species)
FRmm_total <- as.numeric(traindata$FRmm_total)
FRmm_diff <- as.numeric(traindata$FRmm_diff) # Normalized difference
FRmm_3p <-
  as.numeric(paste(traindata$FRmm_3p / FRmm_total)) # Proportion
FRmm_term <-
  as.numeric(paste(traindata$FRmm_term / FRmm_total)) #Proportion
FRmm_AA <-
  as.numeric(paste(traindata$FRmm_AA / FRmm_total)) # Proportion
FRmm_AG <-
  as.numeric(paste(traindata$FRmm_AG / FRmm_total)) # Proportion
FRmm_AC <-
  as.numeric(paste(traindata$FRmm_AC / FRmm_total)) # Proportion
FRmm_TT <-
  as.numeric(paste(traindata$FRmm_TT / FRmm_total)) # Proportion
FRmm_TG <-
  as.numeric(paste(traindata$FRmm_TG / FRmm_total)) # Proportion
FRmm_TC <-
  as.numeric(paste(traindata$FRmm_TC / FRmm_total)) # Proportion
FRmm_GG <-
  as.numeric(paste(traindata$FRmm_GG / FRmm_total)) # Proportion
FRmm_CC <-
  as.numeric(paste(traindata$FRmm_CC / FRmm_total)) # Proportion
FR_length <-
  as.numeric(paste((traindata$F_length + traindata$R_length) / 2)) # Mean
FR_Tm <-
  as.numeric(paste((traindata$F_Tm + traindata$R_Tm) / 2)) # Mean
FR_Tmdiff <-
  as.numeric(paste(abs(traindata$F_Tm - traindata$R_Tm))) # Raw difference
FR_GC <-
  as.numeric(paste((traindata$F_GC + traindata$R_GC) / 2)) # Mean
SYBR_results <- as.factor(traindata$SYBR_results)
SYBR_keep <- as.factor(traindata$SYBR_keep)

traindata <-
  data.frame(
    Assay,
    Species,
    FRmm_total,
    FRmm_diff,
    FRmm_3p,
    FRmm_term,
    FRmm_AA,
    FRmm_AG,
    FRmm_AC,
    FRmm_TT,
    FRmm_TG,
    FRmm_TC,
    FRmm_GG,
    FRmm_CC,
    FR_length,
    FR_Tm,
    FR_Tmdiff,
    FR_GC,
    SYBR_results,
    SYBR_keep
  )

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
traindata[is.nan(traindata)] <- 0

traindata <- subset.data.frame(traindata, SYBR_keep == TRUE)
traindata <- subset.data.frame(traindata, select = c(3:19))

##################################################################################################
### Pre-process SYBR data prior to cross-validation
### Assess response variable distribution
table(traindata$SYBR_results) # Binary outcomes are roughly equal

### Assess highly correlated variables
varcor <- cor(traindata[,-17])
findCorrelation(varcor, cutoff = 0.75) # Length strongly correlated with GC content

### Remove GC content
traindata <-
  subset.data.frame(traindata, select = -16)

### Assess linearly dependent variables
findLinearCombos(traindata[,-16]) # No linearly dependent predictors

### Normalizing data not necessary for random forest models

##################################################################################################
### Implement repeated k-fold cross-validation
traincontrol <-
  trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 10,
    savePredictions = TRUE,
    summaryFunction = prSummary,
    classProbs = TRUE
  )
tunegrid <- expand.grid(.mtry = c(1:15))
set.seed(1111)
trainmodel_sybr <-
  train(
    SYBR_results ~ .,
    data = traindata,
    method = "rf",
    num.trees = 1000,
    metric = "Recall",
    trControl = traincontrol,
    tuneGrid = tunegrid
  )
#save(trainmodel_sybr, file = "SYBR_rf_model.RData")
beep(sound = 1)

##################################################################################################
### Assess model performance
setwd("C:/Users/jkronenberger/Box/ESTCP/eDNAssay/eDNAssay_files_for_GitHub")
load("SYBR_rf_model.RData")
print(trainmodel_sybr)

# trainresults_sybr <- predict(trainmodel_sybr, type="prob") # To output assignment probabilities
trainresults_sybr <- predict(trainmodel_sybr)
traincm_sybr <-
  confusionMatrix(trainresults_sybr, traindata[, 16])
print(traincm_sybr)

varImp(trainmodel_sybr, scale = TRUE)$importance
plot(varImp(trainmodel_sybr, scale = TRUE))

##################################################################################################
### Plot ROC curve
trainmodel_sybr$pred$pred <-
  as.character(trainmodel_sybr$pred$pred)
trainmodel_sybr$pred$obs <-
  as.character(trainmodel_sybr$pred$obs)
trainmodel_sybr$pred <-
  data.frame(trainmodel_sybr$pred, stringsAsFactors = FALSE)
trainmodel_sybr$pred <-
  replace(trainmodel_sybr$pred, trainmodel_sybr$pred == "Amp", 1)
trainmodel_sybr$pred <-
  replace(trainmodel_sybr$pred, trainmodel_sybr$pred == "NoAmp", 0)
index_sybr <- trainmodel_sybr$pred$mtry == 7

trainauc_sybr <-
  ggplot(trainmodel_sybr$pred[index_sybr,], aes(m = Amp, d = as.integer(obs))) +
  geom_roc(n.cuts = 0, color = "black") + coord_equal() + style_roc() +
  ggtitle("ROC") +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) +
  theme(panel.border = element_rect(colour = "black", size = 1.25))

trainauc_sybr <-
  trainauc_sybr + annotate(
    "text",
    x = 0.375,
    y = 0.825,
    size = 5,
    label = paste("AUC =", round((
      calc_auc(trainauc_sybr)
    )$AUC, 3))
  )

print(trainauc_sybr)
