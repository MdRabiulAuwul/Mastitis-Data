
library(DESeq2)

#   Differential expression analysis with limma for GSE15025
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE15025", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL2112", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000001111100000111110000XXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("T","C"))
batch <- pData(GSE15025)$time.point  # or $characteristics_ch1 (check actual metadata)
design <- model.matrix(~ batch + group1 + 0)  
colnames(design) <- gsub("^batch", "B", colnames(design))
colnames(design) <- gsub("^group", "G", colnames(design))
Data_corrected <- removeBatchEffect(scale(Data), batch=batch, design=model.matrix(~ group))
fit <- lmFit(scale(Data), design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="BH", n = nrow(Data), sort.by = "none")
pval <- tT$P.Value
TopDEGT <- which(tT$adj.P.Val < 0.05 & abs(tT$logFC) >= 0.5)

######---------S-S (DESeq2)------##########

colData <- data.frame(condition = group)
ssR <- DESeqDataSetFromMatrix(countData = Data, colData = colData, design = ~ condition)
ssR <- DESeq2::estimateSizeFactors(ssR, type = "ratio")
ssR <- DESeq2::estimateDispersions(ssR, fitType = "parametric")
ssR <- DESeq2::nbinomLRT(ssR, full = ~ condition, reduced = ~ 1)
res <- results(ssR)

#############Machine Learning##########
DataX<-Data[,res] ##Extract only DE gene
library(class)
library(caret)
library(MASS)
library(e1071)
library(rpart)
library(randomForest)
smp_size <- floor(0.75 * nrow(t(DataX)))
## set the seed to make your partition reproducible
train_ind <- sample(seq_len(nrow(t(DataX))), size = smp_size)
DatTra <- DataX[,train_ind]
DatTee <- DataX[,-train_ind]
Trlevel<-group[train_ind]
Televel<-group[-train_ind]

svm.model <- svm(Trlevel ~ ., data = t(DatTra),probability=TRUE)
svm.pred <- predict(svm.model, t(DatTee))

model.nb <- naiveBayes(Trlevel ~ ., data = as.data.frame(t(DatTra)))
pred.nb<-predict(model.nb, as.data.frame(t(DatTee)))

z<-lda(t(DatTra),Trlevel)
p.lda=predict(z, t(DatTee))$class

rm.model <- randomForest(Trlevel ~ ., data =t(DatTra))
predict.rf=predict(rm.model,t(DatTee))

#Repeat the results for 1000 times and averaging the all results
