library(limma)
library(DESeq2)
library(edgeR)

Data1 <- read.csv("GSE64801.csv")
Data<-Data1[,-1]
rownames(Data)<-Data1[,1]

group<-factor(rep(c("1","0"),c(9,5)))
design<-model.matrix(~ group)

######---------S-S (DESeq2)------##########

colData <- data.frame(condition = group)
ssR <- DESeqDataSetFromMatrix(countData = Data, colData = colData, design = ~ condition)
ssR <- DESeq2::estimateSizeFactors(ssR, type = "ratio")
ssR <- DESeq2::estimateDispersions(ssR, fitType = "parametric")
ssR <- DESeq2::nbinomLRT(ssR, full = ~ condition, reduced = ~ 1)
res <- results(ssR)

#############ML##########
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
