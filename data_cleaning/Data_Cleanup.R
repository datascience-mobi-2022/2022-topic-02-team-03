#-------------------------------------------
# general cleaning of datasets: renaming rows, removing duplicates, removing NAs
#-------------------------------------------
library(ggplot2) 
library(biomaRt)

#loading data from data folder + making a copy
tcga_exp <- readRDS("./data/tcga_tumor_log2TPM.RDS")
tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")
tumor_vs_norm <- readRDS("./data/tcga_tumor_normal_datascience_proj_2022.RDS")
genesets <- readRDS("./data/hallmarks_genesets.rds")
tcga_exp_copy <- as.data.frame(t(tcga_exp))

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
par(mar=c(12,2,4,2)+.1)
barplot(NA_sources[-which(NA_sources < 200)], ylim = c(0, 10000), main = "distribution of NA sources in anno data", las=2, sub = "variables with less than 200 NA entries have been omitted", names.arg = )


#set up biomart
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
biotypes <- getBM(attributes = c("gene_biotype","ensembl_gene_id_version"), values = ensemble, mart=ensembl, filters = "ensembl_gene_id_version")

#check biotypes of given genesets: example Prol_AACR

