# (3) Build Prediction models for in vitro data (Dataset = Mad)
# BISC 577A, Unit 3

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(DNAshapeR)
library(caret)
library(e1071)
library(ROCR)

## DNA Shape Prediction
fn_fasta <- "/Users/JC/Desktop/BISC577/gcPBM/Mad.txt.fa"
pred <- getShape(fn_fasta)

## Set feature vectors (1-mer + shape)
featureType6 <- c("1-mer", "1-shape")
featureVector6 <- encodeSeqShape(fn_fasta, pred, featureType6)
head(featureVector6)

## Set feature vectors (1-mer)
featureType5 <- c("1-mer")
featureVector5 <- encodeSeqShape(fn_fasta, pred, featureType5)
head(featureVector5)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/JC/Desktop/BISC577/gcPBM/Mad.txt"
exp_data <- read.table(fn_exp)
df6 <- data.frame(affinity=exp_data$V2, featureVector6)
df5 <- data.frame(affinity=exp_data$V2, featureVector5)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction with L2-regularized (1-mer + shape)
model6 <- train(affinity~., data = df6, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model6

# Average Rsquared (1-mer + shape)
average_Rsquared6 <- mean(na.omit(model6$results$Rsquared))
head(average_Rsquared6)

# Prediction with L2-regularized (1-mer)
model5 <- train(affinity~., data = df5, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model5

# Average Rsquared (1-mer)
average_Rsquared5 <- mean(na.omit(model5$results$Rsquared))
head(average_Rsquared5)

# Rsquared Comparison Plot
plot(model5$results$Rsquared, model6$results$Rsquared, col="blue", xlab="Rsquared (1-mer)", ylab="Rsquared (1-mer+shape)",main ="Comparison of Two Models on Mad")
lines(x = seq(0,1), y = seq(0,1), lty=2)
