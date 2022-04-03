##################################################################################################
### Script trains a random forest classifier to predict SYBR Green-based qPCR cross-amplification
##################################################################################################
library(tidyverse)
library(caret)
library(randomForest)
library(beepr)
library(plotROC)

setwd("C:/Users/jkronenberger/Box/ESTCP/Assays/WORKING")

### Select input CSV file containing training data
traindata_all <-
  read.csv(file.choose(), header = TRUE) # Training data

### Select variables of interest
traindata <-
  subset.data.frame(traindata_all, select = c(4, 5, 7:19, 22))

##################################################################################################
### Pre-process SYBR data prior to cross-validation
### Assess response variable distribution
table(traindata$Results) # Binary outcomes are roughly equal

### Assess highly correlated variables
varcor <- cor(traindata[, -16])
findCorrelation(varcor, cutoff = 0.75) # No highly correlated variables

### Assess linearly dependent variables
findLinearCombos(traindata[, -16]) # No linearly dependent predictors

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
    Results ~ .,
    data = traindata,
    method = "rf",
    num.trees = 1000,
    metric = "Recall",
    trControl = traincontrol,
    tuneGrid = tunegrid
  )
save(trainmodel_sybr, file = "SYBR_trained_model.RData")
beep(sound = 1)

##################################################################################################
### Assess model performance
load("SYBR_trained_model.RData")
print(trainmodel_sybr)

# trainresults_sybr <- predict(trainmodel_sybr, type="prob") # To output assignment probabilities
trainresults_sybr <- predict(trainmodel_sybr)

traindata[, 16] <- as.factor(traindata[, 16])
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
  ggplot(trainmodel_sybr$pred[index_sybr, ], aes(m = Amp, d = as.integer(obs))) +
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
