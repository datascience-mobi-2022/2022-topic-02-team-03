#-------------------------------------------
# general cleaning of datasets: renaming rows, removing unwanted biotypes and low var
#-------------------------------------------
library(ggplot2) 
library(biomaRt)

#loading data from data folder + making a copy
tcga_exp <- readRDS("./data/tcga_tumor_log2TPM.RDS")
tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")
tumor_vs_norm <- readRDS("./data/tcga_tumor_normal_datascience_proj_2022.RDS")
genesets <- readRDS("./data/hallmarks_genesets.rds")
tcga_exp_copy <- as.data.frame(t(tcga_exp))
LUAD_patients <- tcga_exp_copy[which(tcga_anno$cancer_type_abbreviation=="LUAD")]
#----------------------------------------------
#extracting gene IDs and ensemble IDs
IDs <- c()
for (element in rownames(tcga_exp)) {
  IDs <- c(IDs, strsplit(element, "|", fixed = TRUE)[[1]][2])
}

ensemble <- c()
for (element in rownames(tcga_exp)) {
  ensemble <- c(ensemble, strsplit(element, "|", fixed = TRUE)[[1]][1])
}

# setting row names to gene IDs; enumerate duplicates
colnames(tcga_exp_copy) <- make.names(IDs, unique = TRUE)

# check for NAs in both exp and anno
print(paste0("There are " , sum(is.na(tcga_exp)) , " NAs in tcga_exp"))
tcga_anno[tcga_anno == ""] <- NA # some cells were empty but not NA so replace "" with NA
print(paste0("There are " , sum(is.na(tcga_anno)) , " NAs in tcga_anno"))

# where do these NAs come from?
NA_sources <- sort(apply(tcga_anno, 2, function(x){sum(is.na(x))}))
print(paste0("There are " , sum(is.na(tcga_anno$cancer_type_abbreviation)) , " NAs in cancer types"))
par(mar=c(12,4,4,2)+.1)
barplot(NA_sources[-which(NA_sources < 200)], ylim = c(0, 10000), main = "distribution of NA sources in anno data", las=2, names.arg = )

# finding the max variance genes
gene_variances <- sort(apply(tcga_exp_copy, 2, var), decreasing = TRUE)
plot(gene_variances, type="l", xlab = "genes", ylab="variance")

#delete the lowest 40% of genes 
exp_highvar <- tcga_exp_copy[,which(gene_variances > quantile(gene_variances, 0.4))]
gene_variances_highvar <- sort(apply(exp_highvar, 2, var), decreasing = TRUE)
plot(gene_variances_highvar, type="l", xlab = "genes", ylab="variance")

# histogram of mean of genes for LUAD patients 
hist(apply(LUAD_patients, 2, mean),  breaks = 20, xlab = "mean expression over all genes", main="distribution of gene expression in LUAD")

# max and min expressed genes in LUAD
LUAD_means <- apply(LUAD_patients, 2, mean)
LUAD_means_desc <- sort(LUAD_means, decreasing = TRUE)
LUAD_means_asc <- sort(LUAD_means, decreasing = FALSE)
LUAD_means <- data.frame("name" = names(LUAD_means_desc), "desc"=LUAD_means_desc)
rm(LUAD_means_asc, LUAD_means_desc)

# get biotypes with biomart
ensemble_noVersion <- c()
for (element in ensemble) {
  ensemble_noVersion <- c(ensemble_noVersion, strsplit(element, ".", fixed = TRUE)[[1]][1])
}
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biotypes <- getBM(mart = ensembl, values = ensemble_noVersion, filters = "ensembl_gene_id", attributes = c("gene_biotype", "ensembl_gene_id"))

biotypes_paired <- data.frame()
for (gene in biotypes$ensembl_gene_id) {
  biotypes_paired <- rbind(biotypes_paired, tcga_exp[grep(gene, rownames(tcga_exp)),])
  print(grep(toString(gene), rownames(tcga_exp)))
  
}
