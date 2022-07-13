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

##clean and afterwards split dataset into train set and validation set
cor.mat <- cor(t(all_gsva_C5_highSD), method= "spearman")
cor.heatmap <- ComplexHeatmap::Heatmap(cor.mat)

cor.mat_rm <- cor.mat
cor.mat_rm[upper.tri(cor.mat_rm)] <- 0
cor.mat_rm <- as.data.frame(cor.mat_rm)

lowcor_all_gsva <- all_gsva[which(rownames(all_gsva) %in% rownames(differential_genesets[1:50,])),]



# add column to annotation data with classification whether it is LUAD or not
lowcor_all_gsva <- rbind(lowcor_all_gsva, ifelse(tcga_anno$cancer_type_abbreviation == "LUAD", 1, 0))
rownames(lowcor_all_gsva)[51] <- "classifier"
lowcor_all_gsva <- as.data.frame(t(lowcor_all_gsva))

train_dataset <- as.data.frame(lowcor_all_gsva[1:6820,])
test_dataset <- as.data.frame(lowcor_all_gsva[6820:9741,])


#train the model with all pathways
model_all_genesets <- glm(classifier ~  MYELOID_LEUKOCYTE_MIGRATION + REGULATION_OF_AMPA_RECEPTOR_ACTIVITY + T_CELL_MIGRATION + CELLULAR_EXTRAVASATION + POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION + REGULATION_OF_LEUKOCYTE_APOPTOTIC_PROCESS + POSITIVE_REGULATION_OF_INTERLEUKIN_1_BETA_PRODUCTION  + MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE + REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY + REGULATION_OF_B_CELL_PROLIFERATION + B_CELL_DIFFERENTIATION + GRANULOCYTE_ACTIVATION + NEGATIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION + REGULATION_OF_PHAGOCYTOSIS + POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION + POSITIVE_REGULATION_OF_LEUKOCYTE_CHEMOTAXIS + REGULATION_OF_LEUKOCYTE_MIGRATION + REGULATION_OF_LEUKOCYTE_CHEMOTAXIS, family = binomial(link = "logit"), data = train_dataset)
blr_step_aic_both(model_all_genesets)
plot(blr_step_aic_both(model_all_genesets))

blr_regress(model_all_genesets)

# check for fitting cutoff value
predict_train <- predict(model_all_genesets, type = "response")
tapply(predict_train, train_dataset$classifier, mean)

# test
predict_test <- predict(model_all_genesets, type = "response", newdata = test_dataset)
addmargins(table(test_dataset$classifier, predict_test > 0.35))

#model evaluation
blr_model_fit_stats(model_all_genesets)
