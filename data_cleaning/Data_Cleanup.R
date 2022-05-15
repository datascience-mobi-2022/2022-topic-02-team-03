#-------------------------------------------
# general cleaning of datasets: renaming rows, removing unwanted biotypes and low var
#-------------------------------------------
library(ggplot2) 
library(biomaRt)
library(parallel)
library(grid)


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

rm(element)
# setting row names to gene IDs; enumerate duplicates
colnames(tcga_exp_copy) <- make.names(IDs, unique = TRUE)

# check for NAs in both exp and anno
print(paste0("There are " , sum(is.na(tcga_exp)) , " NAs in tcga_exp"))
tcga_anno[tcga_anno == ""] <- NA # some cells were empty but not NA so replace "" with NA
print(paste0("There are " , sum(is.na(tcga_anno)) , " NAs in tcga_anno"))

# where do these NAs come from?
NA_sources <- sort(apply(tcga_anno, 2, function(x){sum(is.na(x))}))
NA_sources <- data.frame("sum" = NA_sources, "variable" = names(NA_sources))
print(paste0("There are " , sum(is.na(tcga_anno$cancer_type_abbreviation)) , " NAs in cancer types"))

ggplot(NA_sources, aes(x=reorder(variable, sum), y = sum))+
  geom_col(color = "dodgerblue3")+
  coord_flip()+
  xlab("")+
  ylab("sum of NAs in anno")


# max and min expressed genes in LUAD
LUAD_means <- apply(LUAD_patients, 2, mean)
LUAD_means_desc <- sort(LUAD_means, decreasing = TRUE)
LUAD_means_asc <- sort(LUAD_means, decreasing = FALSE)
LUAD_means <- data.frame("name" = names(LUAD_means_desc), "desc"=LUAD_means_desc)
all_means <- apply(tcga_exp_copy, 2, mean)
all_means <- data.frame("name" = names(all_means), "desc" = sort(all_means, decreasing = TRUE))
# histogram of mean of genes for LUAD patients 
#hist(apply(LUAD_patients, 2, mean),  breaks = 20, xlab = "mean expression over all genes", main="distribution of gene expression in LUAD")

ggplot(all_means, aes(x = desc))+
  geom_histogram(bins=30, color = "grey")+
  scale_y_log10()+
  xlab("mean expression of gene")+
  ylab("log10 count")+
  ggtitle("distribution of mean gene expression in all patients")

rm(LUAD_means_asc, LUAD_means_desc)
# get biotypes of exp-genes with biomart
ensemble_noVersion <- c()
for (element in ensemble) {
  ensemble_noVersion <- c(ensemble_noVersion, strsplit(element, ".", fixed = TRUE)[[1]][1])
}
rm(element)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biotypes <- getBM(mart = ensembl, values = ensemble_noVersion, filters = "ensembl_gene_id", attributes = c("gene_biotype", "ensembl_gene_id"))

genes_with_biotypes <- list()
genes_with_biotypes <- mclapply(biotypes$ensembl_gene_id, function(gene){grep(gene, rownames(tcga_exp))}, mc.cores = detectCores())
genes_with_biotypes <- unlist(genes_with_biotypes)

genes_and_biotypes <- data.frame("ensembl_ID" = ensemble[genes_with_biotypes], "gene_name" = colnames(tcga_exp_copy)[genes_with_biotypes], "biotype" = biotypes$gene_biotype)

# looking up biotypes of the hallmarks to delete all genes that do not fit biotype
biotypes_hallmarks <- getBM(mart = ensembl, values = unique(unlist(genesets$genesets, use.names = FALSE)), filters = "hgnc_symbol", attributes = c("gene_biotype", "ensembl_gene_id", "hgnc_symbol"))
relevant_biotypes <- levels(as.factor(biotypes_hallmarks$gene_biotype))

biotypes_metabol <- getBM(mart = ensembl, values = unique(unlist(metabolism_list_gs, use.names = FALSE)), filters = "hgnc_symbol", attributes = c("gene_biotype", "ensembl_gene_id", "hgnc_symbol"))
relevant_biotypes <- unique(append(relevant_biotypes, levels(as.factor(biotypes_metabol$gene_biotype))))

# delete all genes that do not fit biotype in biotypes df
biotypes <- biotypes[which(biotypes$gene_biotype %in% relevant_biotypes),]

#make tcga_copy smaller by deleting wrong biotype genes
tcga_exp_copy <- tcga_exp_copy[,which(ensemble_noVersion %in% biotypes$ensembl_gene_id)]
IDs <- c()
for (element in colnames(tcga_exp_copy)) {
  IDs <- c(IDs, strsplit(element, "|", fixed = TRUE)[[1]][2])
}
rm(element)
colnames(tcga_exp_copy) <- make.names(IDs, unique = TRUE)

# finding the max variance genes
gene_variances <- sort(apply(tcga_exp_copy, 2, var), decreasing = TRUE)
#plot(gene_variances, type="l", xlab = "genes", ylab="variance")

#delete the lowest 50% of genes 
exp_highvar <- tcga_exp_copy[,which(gene_variances > quantile(gene_variances, 0.5))]
gene_variances_highvar <- sort(apply(exp_highvar, 2, var), decreasing = TRUE)


grob <- grobTree(textGrob("AL162151.3", x=0.88,  y=0.09, hjust=0, gp=gpar(col="red", fontsize=10, fontface="italic")))

ggplot(as.data.frame(gene_variances_highvar), aes(gene_variances_highvar))+
  geom_freqpoly(bins=50) +
  scale_y_log10() +
  ggtitle("Distribution of variances in TCGA_exp")+
  xlab("variance of gene")+
  annotation_custom(grob)
  

saveRDS(exp_highvar, "./data/tcga_exp_small.RDS")


