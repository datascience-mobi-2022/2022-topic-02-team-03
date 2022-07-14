#logistic regression on tcga_exp to evaluate tumor type form pathway expression data
## load needed packages
library(rms)
library(MASS)
library(DAAG)
library(pROC)
library(ComplexHeatmap)
library(caret)
library(Seurat)
library(blorr)
library(ggplot2)
library(ROCR)

tcga_exp_short <- readRDS("./data/tcga_exp_small.RDS")
tcga_anno <- readRDS(".data/tcga_tumor_annotation.RDS")
##clean and afterwards split dataset into train set and validation set
# cor.mat <- cor(t(all_gsva_C5_highSD), method= "spearman")
# cor.heatmap <- ComplexHeatmap::Heatmap(cor.mat)
# 
# cor.mat_rm <- cor.mat
# cor.mat_rm[upper.tri(cor.mat_rm)] <- 0
# cor.mat_rm <- as.data.frame(cor.mat_rm)
# 
# lowcor_all_gsva <- all_gsva[which(rownames(all_gsva) %in% rownames(differential_genesets[1:50,])),]

# extract most diverse gene between LUAD and the other cancer types
tcga_luad <- tcga_exp_short[,which(tcga_anno$cancer_type_abbreviation == "LUAD")]
tcga_rest <- tcga_exp_short[,-which(tcga_anno$cancer_type_abbreviation == "LUAD")]

luad_means <- apply(tcga_luad, 1, mean)
rest_means <- apply(tcga_rest, 1, mean)

luad_rest_diff <- luad_means - rest_means

gene_candidates <- c(sort(luad_rest_diff, decreasing = TRUE)[1:10], sort(luad_rest_diff, decreasing = FALSE)[1:10])
gene_candidates <- names(gene_candidates)

#qc run PCAUMAP on 20 candidates
QC_pcaumap <- as.data.frame(umap(RunPCA(as.matrix(tcga_exp_short[gene_candidates,]), npcs = 5)@cell.embeddings, metric = "cosine"))

ggplot(QC_pcaumap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()



# add column to annotation data with classification whether it is LUAD or not
lowcor_all_gsva <- rbind(tcga_exp_short, ifelse(tcga_anno$cancer_type_abbreviation == "LUAD", 1, 0))
rownames(lowcor_all_gsva)[15650] <- "classifier"
lowcor_all_gsva <- as.data.frame(t(lowcor_all_gsva))

train_dataset <- as.data.frame(lowcor_all_gsva[1:6820,])
test_dataset <- as.data.frame(lowcor_all_gsva[6820:9741,])


#train the model with genes
model_all_genesets <- glm(classifier ~   SFTPB + EN1 + BPIFA1 + C10orf99 + C4BPA + HAND2, family = binomial(link = "logit"), data = train_dataset)
blr_step_aic_both(model_all_genesets)
plot(blr_step_aic_both(model_all_genesets))

blr_regress(model_all_genesets)

# check for fitting cutoff value
predict_train <- predict(model_all_genesets, type = "response")
tapply(predict_train, train_dataset$classifier, mean)

# test
predict_test <- predict(model_all_genesets, type = "response", newdata = test_dataset)
predict_test <- ifelse(predict_test > 0.5, 1, 0)
addmargins(table(test_dataset$classifier, predict_test))

#model evaluation
blr_model_fit_stats(model_all_genesets)

missing_classerr <- mean(predict_test != test_dataset$classifier)
print(paste('Accuracy =', 1 - missing_classerr))

ROCPred <- prediction(predict_test, test_dataset$classifier)
ROCPer <- performance(ROCPred, measure = "tpr", x.measure = "fpr")


auc <- performance(ROCPred, measure = "auc")
auc <- auc@y.values[[1]]
auc

plot(ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC CURVE")
abline(a = 0, b = 1)
auc <- round(auc, 4)
legend(.6, .4, auc, title = "AUC", cex = 1)
