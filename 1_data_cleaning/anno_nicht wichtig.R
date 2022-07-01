#----------------------------------------------
#loading data
tcga_exp <- readRDS("./data/tcga_tumor_log2TPM.RDS")
tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")
tumor_vs_norm <- readRDS("./data/tcga_tumor_normal_datascience_proj_2022.RDS")
genesets <- readRDS("./data/hallmarks_genesets.rds")




#----------------------tcga_anno------------------------
# create a new dataframe with only LUAD Patient Info
tcga_anno_split = split(tcga_anno, tcga_anno$cancer_type_abbreviation)          #splits the data frame into categories of can._typ._abbr.

tcga_anno_split_luad = tcga_anno_split$LUAD


# df only with LUAD patients
 ggplot(anno_LUAD_patients, aes(age_at_initial_pathologic_diagnosis, gender, col = gender))+
   geom_point()


anno_LUAD_BRCA = rbind(tcga_anno_split$BRCA, tcga_anno_split$LUAD)# df LUAD and BRCA patients


#--------genesets-------------
class(genesets)
names(genesets)
names(genesets$genesets)

pw = which(names(genesets$genesets) == "Prol_MTOR")# which pathway of the list has the name "..."
genesets$description[pw]# reading the description of this pathway
