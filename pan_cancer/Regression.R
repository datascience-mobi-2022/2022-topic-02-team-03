#logistic regression on tcga_exp to evaluate tumor type form pathway expression data
## load needed packages
library(rms)
library(MASS)
library(DAAG)
library(pROC)
library(ComplexHeatmap)
library(caret)

##clean and afterwards split dataset into train set and validation set
cor.mat <- cor(t(all_gsva), method= "spearman")
cor.heatmap <- ComplexHeatmap::Heatmap(cor.mat)

cor.mat_rm <- cor.mat
cor.mat_rm[upper.tri(cor.mat_rm)] <- 0
cor.mat_rm <- as.data.frame(cor.mat_rm)

lowcor_all_gsva <- as.data.frame(all_gsva[-findCorrelation(cor.mat),])



train_dataset <- lowcor_all_gsva[,1:6819]
test_dataset <- lowcor_all_gsva[,6820:9741]

