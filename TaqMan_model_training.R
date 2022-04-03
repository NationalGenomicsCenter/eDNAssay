##################################################################################################
### Script trains a random forest classifier to predict TaqMan-based qPCR cross-amplification
##################################################################################################
library(tidyverse)
library(smotefamily)
library(caret)
library(randomForest)
library(beepr)
library(plotROC)

### Select input CSV file containing training data
traindata_all <-
  read.csv(file.choose(), header = TRUE) # Training data

### Select variables of interest
traindata <-
  subset.data.frame(traindata_all, select = c(4, 5, 7:19, 21:32, 34))

##################################################################################################
### Pre-process TaqMan data prior to cross-validation
### Assess response variable distribution
table(traindata$Results) # Binary outcomes not very equal

### Balance class variables using SMOTE
traindatas <-
  SMOTE(traindata[, -28], traindata[, 28], dup_size = 2)
traindatas <- traindatas$data
traindatas$class <- as.factor(traindatas$class)
table(traindatas$class)

### Assess variables with zero or near zero-variance (can affect model performance, but random forest impervious to them)
nearZeroVar(traindatas, saveMetrics = TRUE)

### Assess highly correlated variables
varcor <- cor(traindatas[,-28])
findCorrelation(varcor, cutoff = 0.75) # P_Tm correlated with F_Tmprop but okay

### Assess linearly dependent variables
findLinearCombos(traindatas[, -28]) # No linearly dependent variables

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
tunegrid <- expand.grid(.mtry = c(1:27))
set.seed(1111)
trainmodel_taqman <-
  train(
    class ~ .,
    data = traindatas,
    method = "rf",
    ntree = 1000,
    metric = "Recall",
    trControl = traincontrol,
    tuneGrid = tunegrid
  )
save(trainmodel_taqman, file = "TaqMan_trained_model.RData")
beep(sound = 1)

##################################################################################################
### Assess model performance
load("TaqMan_trained_model.RData")
print(trainmodel_taqman)

# trainresults_taqman <- predict(trainmodel_taqman, type="prob") # To output assignment probabilities
trainresults_taqman <- predict(trainmodel_taqman)
traincm_taqman <-
  
  confusionMatrix(trainresults_taqman, traindatas[, 28])
print(traincm_taqman)

varImp(trainmodel_taqman, scale = TRUE)$importance
plot(varImp(trainmodel_taqman, scale = TRUE))

##################################################################################################
### Plot ROC curve
trainmodel_taqman$pred$pred <-
  as.character(trainmodel_taqman$pred$pred)
trainmodel_taqman$pred$obs <-
  as.character(trainmodel_taqman$pred$obs)
trainmodel_taqman$pred <-
  data.frame(trainmodel_taqman$pred, stringsAsFactors = FALSE)
trainmodel_taqman$pred <-
  replace(trainmodel_taqman$pred, trainmodel_taqman$pred == "Amp", 1)
trainmodel_taqman$pred <-
  replace(trainmodel_taqman$pred, trainmodel_taqman$pred == "NoAmp", 0)
index_taqman <- trainmodel_taqman$pred$mtry == 2

trainauc_taqman <-
  ggplot(trainmodel_taqman$pred[index_taqman,], aes(m = Amp, d = as.integer(obs))) +
  geom_roc(n.cuts = 0, color = "black") + coord_equal() + style_roc() +
  ggtitle("ROC") +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) +
  theme(panel.border = element_rect(colour = "black", size = 1.25))

trainauc_taqman <-
  trainauc_taqman + annotate(
    "text",
    x = 0.375,
    y = 0.825,
    size = 5,
    label = paste("AUC =", round((
      calc_auc(trainauc_taqman)
    )$AUC, 3))
  )

print(trainauc_taqman)
