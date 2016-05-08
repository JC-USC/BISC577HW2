# (3) Build Prediction models for in vitro data (Dataset = Max)
# BISC 577A, Unit 3

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(DNAshapeR)
library(caret)
library(e1071)
library(ROCR)

## DNA Shape Prediction
fn_fasta <- "/Users/JC/Desktop/BISC577/gcPBM/Max.txt.fa"
pred <- getShape(fn_fasta)

## Set feature vectors (1-mer + shape)
featureType4 <- c("1-mer", "1-shape")
featureVector4 <- encodeSeqShape(fn_fasta, pred, featureType4)
head(featureVector4)

## Set feature vectors (1-mer)
featureType3 <- c("1-mer")
featureVector3 <- encodeSeqShape(fn_fasta, pred, featureType3)
head(featureVector3)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/JC/Desktop/BISC577/gcPBM/Max.txt"
exp_data <- read.table(fn_exp)
df4 <- data.frame(affinity=exp_data$V2, featureVector4)
df3 <- data.frame(affinity=exp_data$V2, featureVector3)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction with L2-regularized (1-mer + shape)
model4 <- train(affinity~., data = df4, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model4

# Average Rsquared (1-mer + shape)
average_Rsquared4 <- mean(na.omit(model4$results$Rsquared))
head(average_Rsquared4)

# Prediction with L2-regularized (1-mer)
model3 <- train(affinity~., data = df3, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model3

# Average Rsquared (1-mer)
average_Rsquared3 <- mean(na.omit(model3$results$Rsquared))
head(average_Rsquared3)

# Rsquared Comparison Plot
plot(model3$results$Rsquared, model4$results$Rsquared, col="blue", xlab="Rsquared (1-mer)", ylab="Rsquared (1-mer+shape)",main ="Comparison of Two Models on Max")
lines(x = seq(0,1), y = seq(0,1), lty=2)

# Rsquared Comparison Plot (Mad, Max, Myc)
plot(model3$results$Rsquared, model4$results$Rsquared, col="blue", xlab="Rsquared (1-mer)", ylab="Rsquared (1-mer+shape)",main ="Comparison of Two Models on Mad, Max, Myc")
points(model1$results$Rsquared, model2$results$Rsquared, col="red")
points(model5$results$Rsquared, model6$results$Rsquared, col="green")
lines(x = seq(0,1), y = seq(0,1), lty=2)
legend(0.525, 0.8, c("Mad","Max","Myc"), lwd=5, cex=1, col=c("green", "blue", "red"), horiz = TRUE)



