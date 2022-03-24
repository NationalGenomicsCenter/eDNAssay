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
traindata <- read.csv(file.choose(), header = TRUE) # Training data

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
Pmm_total <- as.numeric(traindata$Pmm_total)
Pmm_center <-
  as.numeric(paste(traindata$Pmm_center / Pmm_total)) # Proportion
Pmm_AA <-
  as.numeric(paste(traindata$Pmm_AA / Pmm_total)) # Proportion
Pmm_AG <-
  as.numeric(paste(traindata$Pmm_AG / Pmm_total)) # Proportion
Pmm_AC <-
  as.numeric(paste(traindata$Pmm_AC / Pmm_total)) # Proportion
Pmm_TT <-
  as.numeric(paste(traindata$Pmm_TT / Pmm_total)) # Proportion
Pmm_TG <-
  as.numeric(paste(traindata$Pmm_TG / Pmm_total)) # Proportion
Pmm_TC <-
  as.numeric(paste(traindata$Pmm_TC / Pmm_total)) # Proportion
Pmm_GG <-
  as.numeric(paste(traindata$Pmm_GG / Pmm_total)) # Proportion
Pmm_CC <-
  as.numeric(paste(traindata$Pmm_CC / Pmm_total)) # Proportion
P_length <- as.numeric(traindata$P_length)
P_Tm <- as.numeric(traindata$P_Tm)
P_GC <- as.numeric(traindata$P_GC)
TaqMan_results <- as.factor(traindata$TaqMan_results)

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
    Pmm_total,
    Pmm_center,
    Pmm_AA,
    Pmm_AG,
    Pmm_AC,
    Pmm_TT,
    Pmm_TG,
    Pmm_TC,
    Pmm_GG,
    Pmm_CC,
    P_length,
    P_Tm,
    P_GC,
    TaqMan_results
  )

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
traindata[is.nan(traindata)] <- 0

traindata <-
  subset.data.frame(traindata, select = c(3:32))

##################################################################################################
### Pre-process TaqMan data prior to cross-validation
### Assess response variable distribution
table(traindata$TaqMan_results) # Binary outcomes not very equal

### Balance class variables using SMOTE
traindata_taqmans <-
  SMOTE(traindata[, -30], traindata[, 30], dup_size = 2)
traindata_taqmans <- traindata_taqmans$data
traindata_taqmans$class <- as.factor(traindata_taqmans$class)
table(traindata_taqmans$class)

### Assess variables with zero or near zero-variance (can affect model performance, but random forest impervious to them)
nearZeroVar(traindata_taqmans, saveMetrics = TRUE)

### Assess highly correlated variables
varcor <- cor(traindata_taqmans[,-30])
findCorrelation(varcor, cutoff = 0.75) # Length strongly correlated with GC content

### Remove GC content
traindata_taqmans <-
  subset.data.frame(traindata_taqmans, select = -c(16, 29))

### Assess linearly dependent variables
findLinearCombos(traindata_taqmans[, -28]) # No linearly dependent variables

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
    data = traindata_taqmans,
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
  
  confusionMatrix(trainresults_taqman, traindata_taqmans[, 28])
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
