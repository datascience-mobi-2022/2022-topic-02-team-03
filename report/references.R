# load all packages that were used over the course of the project
library(bayesbio)
library(babyplots)
library(biomaRt)
library(cinaR)
library(cluster)
library(edgeR)
library(enrichplot)
library(FactoMineR)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtools)
library(limma)
library(knitr)
library(msigdbr)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(Seurat)
library(tidyverse)
library(uwot)

#create a .bib file of the citation of the aforementioned packages

packages.bib = write_bib( c (x = .packages(), file = ''))
packages.bib

