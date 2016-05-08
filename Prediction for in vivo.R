## BISC 577, Unit 3
## Prediction for in vivo

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(DNAshapeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Mmusculus.UCSC.mm9)
library(caret)
library(e1071)
library(ROCR)

seqLength <- 30
sampleSize <- 1000

# Bound (ChIP-seq)
ah <- AnnotationHub()
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = ctcfPeaks, BSgenome.Mmusculus.UCSC.mm10, width = seqLength, filename = "bound.fa" )

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()
while ( length(randomGr) < sampleSize ) {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:(tmpChrLength-seqLength),1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}

# Check overlap
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, filename = "unbound.fa")


## Merge bound and unbound data
# Combine two datasets and generate one file
boundFasta <- readBStringSet("bound.fa")
boundFasta <- sample(boundFasta, sampleSize) # Only randomly choose fixed size of data for sampling
unboundFasta <- readBStringSet("unbound.fa")
names(unboundFasta) <- paste0( names(unboundFasta), "_unbound")
writeXStringSet( c(boundFasta, unboundFasta), "ctcf.fa" )

# Generate binding classfication file
boundTxt <- cbind( sapply(1:length(boundFasta), 
                          function(x) as.character(boundFasta[[x]])), 
                   matrix(1, length(boundFasta), 1))
unboundTxt <- cbind( sapply(1:length(unboundFasta),
                            function(x) as.character(unboundFasta[[x]])),
                     matrix(0, length(unboundFasta), 1))
write.table(rbind(boundTxt, unboundTxt), "ctcf.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

## DNAshapeR prediction
shapePred <- getShape("ctcf.fa")

## Set features (1-mer & 1-shape)
featureType1m <- c("1-mer")
featureVector1m <- encodeSeqShape("ctcf.fa", shapePred, featureType1m)
featureType1ms <- c("1-mer", "1-shape")
featureVector1ms <- encodeSeqShape("ctcf.fa", shapePred, featureType1ms)

## Perform logistic regression
exp_data <- read.table("ctcf.txt")

# Prepare data
exp_data$V2 <- ifelse(exp_data$V2 == 1 , "Y", "N")
df1m <- data.frame(isBound = exp_data$V2, featureVector1m)
df1ms <- data.frame(isBound = exp_data$V2, featureVector1ms)

# 2-fold cross-validation (set Caret parameters)
trainControl <- trainControl(method = "cv", number = 2, savePredictions = TRUE, classProbs = TRUE)
model1m <- train(isBound~ ., data = df1m, trControl = trainControl, method = "glm", family = binomial, metric ="ROC")
summary(model1m)
model1ms <- train(isBound~ ., data = df1ms, trControl = trainControl, method = "glm", family = binomial, metric ="ROC")
summary(model1ms)

## Plot AUROC
prediction1m <- prediction( model1m$pred$Y, model1m$pred$obs )
performance1m <- performance( prediction1m, "tpr", "fpr" )
plot(performance1m)
prediction1ms <- prediction( model1ms$pred$Y, model1ms$pred$obs )
performance1ms <- performance( prediction1ms, "tpr", "fpr" )
plot(performance1ms)

auc1m <- performance(prediction1m, "auc")
auc1m <- unlist(slot(auc1m, "y.values"))
auc1m

auc1ms <- performance(prediction1ms, "auc")
auc1ms <- unlist(slot(auc1ms, "y.values"))
auc1ms
