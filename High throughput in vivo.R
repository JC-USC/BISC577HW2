# BISC 577, Unit 3
# High-throughput in vivo

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(DNAshapeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(caret)
library(e1071)
library(ROCR)

## Get Chip-seq data, store in fasta
ah <- AnnotationHub()
ah
unique(ah$dataprovider)
unique(ah$species)

mmus <- ah[["AH28451"]]
seqlevelsStyle(mmus) <- "UCSC"

getFasta(mmus, BSgenome.Mmusculus.UCSC.mm10, width = 400, filename = "mmus.fa")


# Predict DNA shapes
fn <- "mmus.fa"
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$HelT)

heatShape(pred$MGW, 20)
heatShape(pred$ProT, 20)
heatShape(pred$HelT, 20)
