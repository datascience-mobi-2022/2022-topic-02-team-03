#logistic regression on tcga_exp to evaluate tumor type form pathway expression data
## load needed packages
library(rms)
library(MASS)
library(DAAG)
library(pROC)
library(ComlpexHeatmap)

##clean and afterwards split dataset into train set and validation set
cor.mat <- cor(t(all_gsva), method= "spearman")



