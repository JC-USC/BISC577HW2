# (3) Build Prediction models for in vitro data (Dataset = Myc)
# BISC 577A, Unit 3

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(DNAshapeR)
library(caret)
library(e1071)
library(ROCR)

## DNA Shape Prediction
fn_fasta <- "/Users/JC/Desktop/BISC577/gcPBM/Myc.txt.fa"
pred <- getShape(fn_fasta)

## Set feature vectors (1-mer + shape)
featureType2 <- c("1-mer", "1-shape")
featureVector2 <- encodeSeqShape(fn_fasta, pred, featureType2)
head(featureVector2)

## Set feature vectors (1-mer)
featureType1 <- c("1-mer")
featureVector1 <- encodeSeqShape(fn_fasta, pred, featureType1)
head(featureVector1)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/JC/Desktop/BISC577/gcPBM/Myc.txt"
exp_data <- read.table(fn_exp)
df2 <- data.frame(affinity=exp_data$V2, featureVector2)
df1 <- data.frame(affinity=exp_data$V2, featureVector1)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction with L2-regularized (1-mer + shape)
model2 <- train(affinity~., data = df2, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2

# Average Rsquared (1-mer + shape)
average_Rsquared2 <- mean(na.omit(model2$results$Rsquared))
head(average_Rsquared2)

# Prediction with L2-regularized (1-mer)
model1 <- train(affinity~., data = df1, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model1

# Average Rsquared (1-mer)
average_Rsquared1 <- mean(na.omit(model1$results$Rsquared))
head(average_Rsquared1)

# Rsquared Comparison Plot
plot(model1$results$Rsquared, model2$results$Rsquared, col="blue", xlab="Rsquared (1-mer)", ylab="Rsquared (1-mer+shape)",main ="Comparison of Two Models on Myc")
lines(x = seq(0,1), y = seq(0,1), lty=2)
