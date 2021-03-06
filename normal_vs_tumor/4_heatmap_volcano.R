
library(ComplexHeatmap)
library(pheatmap)
library(GSVA)
library(ggplot2)
library("gridExtra")
library(gtools)
library(parallel)
library(grid)
library(ggrepel)
library(EnhancedVolcano)

#GSVA c5 results -> volcano & heatmap

all_gsva_tvn <- readRDS("./data/all_gsva_tvn.RDS")
sd = c(apply(all_gsva_tvn, 1, sd))
q = quantile(sd, probs = c(.8,.9,.95,.98,1), names = F)
q



all_gsva_tvn_cut = as.data.frame(all_gsva_tvn[which(sd >= q[3]),])



#overlapping pathways
#all_gsva_tvn_cut = all_gsva_tvn_cut[na.omit((match(names(genesets$genesets),rownames(all_gsva_tvn_cut)))),]



#renaming the patients to a easier format
t_name = c(1:58)
n_name = c(1:58)

for(i in 1:58){
  t_name[i] = c(paste("t",i))
  n_name[i] = c(paste("n",i))
}

colnames(all_gsva_tvn_cut)[1:58] = n_name
colnames(all_gsva_tvn_cut)[59:116] = t_name

dim(all_gsva_tvn_cut)

Heatmap(all_gsva_tvn_cut,
        row_names_gp =  grid::gpar(fontsize=7),
        column_names_gp = grid::gpar(fontsize = 7),
        row_names_side = "left",
        show_row_dend = F,
        heatmap_legend_param = list(title = "relative expression"),
        row_names_max_width = unit(14, "cm"),
        show_column_names = T
        )
        

# volcano plots
## 1. calculate p values
pw = dim(all_gsva_tvn_cut)[1]

all_gsva_n = as.data.frame(t(all_gsva_tvn_cut[,1:58]))
all_gsva_t = as.data.frame(t(all_gsva_tvn_cut[,59:116]))

p.wilc = c(1:pw)
for(i in 1:pw){
  p.wilc[i] =wilcox.test(all_gsva_n[,i],all_gsva_t[,i], paired = TRUE, alternative = "two.sided", exact = F, correct = F )$p.value
}

mean_n = colMeans(all_gsva_n)

mean_t = colMeans(all_gsva_t)

fc = (mean_t-mean_n)


volc_df = data.frame(p.wilc = p.wilc,
                     mean_n = mean_n,
                     mean_t = mean_t,
                     fc = (fc),
                     names = rownames(all_gsva_tvn_cut))

rownames(volc_df) = rownames(all_gsva_tvn_cut)

#sign p values
for (i in 1:nrow(volc_df)){ #setting up expression
  ifelse((   (volc_df$fc[i] > 0)   ),
         volc_df$spval[i] <- -log10(volc_df$p.wilc[i]),
         volc_df$spval[i] <- log10(volc_df$p.wilc[i]))
}
rm(i)

#sort df by spval
volc_df = volc_df[rev(order(volc_df$fc)),]

alpha = 0.5
alpha.corr = -log10(alpha/nrow(volc_df))
alpha.corr



for (i in 1:nrow(volc_df)){ #setting up expression
  ifelse((   (volc_df$fc[i] > 0) & (-log10(volc_df$p.wilc[i]) > alpha.corr)   ), volc_df$expr[i] <- "up",
         ifelse((   (volc_df$fc[i] < 0) & (-log10(volc_df$p.wilc[i]) > alpha.corr)   ),volc_df$expr[i] <- "down",
                volc_df$expr[i] <- "not sign."))
}

selection1 = volc_df[c(1,4,7,13,21,35,116,114,84,66, 62),]
selection1 = volc_df[which(volc_df$expr == "down"),]
selection1 = volc_df[which(volc_df$expr == "not sign."),]
selection1 = volc_df[which(volc_df$expr == "up"),]
selection1 = rbind(head(volc_df, 5), tail(volc_df, 3))
selection1 = volc_df[which(volc_df$names == "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"),]





volcano = ggplot(volc_df, aes(fc, -log10(p.wilc), col = expr))+ 
  geom_point(size = 4, alpha = .4)+
  labs(title = "pathway expression", y = "-log10(p.value)", x = "log2(Foldchange)")+
  geom_hline(yintercept = alpha.corr,
             linetype = "dashed")+
  scale_color_manual(values=c("blue", "gray48", "red"))+
  geom_point(data = selection1,
             aes(fc, -log10(p.wilc)), size = 2, col = "green", alpha = .9)+
  geom_label_repel(data = selection1,
                  aes(fc, -log10(p.wilc), label = rownames(selection1)),
                  size = 3, hjust = -.1, col = "black", fill = "white", alpha = .8,
                  min.segment.length = 0
                  )


volcano + theme_light()



#costumized volcano

#select genesets
all_final_genesets_C5 <- readRDS("./data/all_final_genesets(C5).RDS")


genesets = all_final_genesets_C5[match(rownames(volc_df), names(all_final_genesets_C5))]


LUAD_tumor_vs_norm_clean <- readRDS("./data/LUAD_tumor_vs_norm_clean.RDS")
luad.info = LUAD_tumor_vs_norm_clean[["infomatrix"]]
rm(LUAD_tumor_vs_norm_clean)



















alpha = 0.05
alpha.corr = -log10(alpha/nrow(luad.info))


z = length(genesets)
pw.v.plots = list()

name = c()
regulation = data.frame(up = c(1:z),
                        down = c(1:z),
                        not_sign = c(1:z),
                        number_of_genes = c(1:z))


which(names(genesets) == "CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS")
for(i in 1:z){
  # select genes from one pwathway
  x = 116 # choose pathway by position
  
  pw = names(genesets[x])
  pw
  
  pw.genes = unique(unlist(genesets[x], use.names = F))#genes in the pathway that are at the leading edge
  
  genes = luad.info$name #all genes
  positions = sort(match(pw.genes, genes))#compare to get positions
  positions
  
  selection = luad.info[positions,] # select rows from data frame
  selection$pathway = pw # make new column, called like the pw
  
  sel = luad.info$diff.mean[positions]#select mean expr for genes in pw
  mean = mean(sel)  # calculate mean expression of the pathway the genes used are the ones in the leading edge
  median = median(sel)
  
  
  
  
  #calculate % genes by up or down regulation
  hits = paste(dim(selection)[1], "genes")
  p.down = (length(which(selection$expr == "down"))/dim(selection)[1])*100
  p.up = (length(which(selection$expr == "up"))/dim(selection)[1])*100
  
  hits
  percent.down = paste(round(p.down, digits = 2), "% down")
  percent.up = paste(round(p.up, digits = 2), "% up")
  rest = 100 - round(p.down, digits = 2) - round(p.up, digits = 2)
  ns = paste(rest, "% not sign.")
  vh = round(p.up/p.down, digits = 2)
  
  #calculate for each pw the menan expression of the genes
  
  
  #
  ##
  ####
  ##
  #
  selection2 = selection[(order(selection$diff.med)),]
  selection2 = head(selection2,5)
  
  col = c("gray0", "firebrick3", "deepskyblue", "gray75") # choosing the colours
  ###plots
  #1 ) normaler v.plot
  plot.v <- ggplot(luad.info, aes(y = abs(spval), diff.mean, colour = expr, label = name))+
    geom_point(size = 2, alpha = .4)+
    labs(title = paste(pw), x = "log2(fc)", y = "-log10(pwilc)")+
    scale_color_manual(name = "",
                       breaks = c(pw, "up", "down", "not sign."),
                       labels = c(hits, percent.up, percent.down, ns),
                       values = col )+ 
    geom_point(data = selection,
               aes(y = abs(spval), diff.mean, color = pathway), size = 3)+
    geom_label_repel(data = selection2, aes(y = abs(spval), diff.mean, color = pathway), size = 5,
                     col = "black",
                     fill = "white", alpha = .8,
                     min.segment.length = 0)+
    theme_light()
  
  
  ###saving
  pw.v.plots[[i]] <- plot.v #save the plot in the list element

  message("plot saved")
  
  regulation$up[i] = round(p.up, digits = 2)
  regulation$down[i] = round(p.down, digits = 2)
  regulation$not_sign[i] = rest
  regulation$up_div_down[i] = vh
  regulation$number_of_genes[i] = dim(selection)[1]
  regulation$mean[i] = mean
  regulation$median[i] = median
  
  message("regulation filled")
  
  name[i] <- paste(pw) #names for list elements
  message("name saved",i)
  
  
  
}
names(pw.v.plots) <- paste0(name)

#new list object with all plots
plots = list(pw.v.plots = pw.v.plots,
             regulation = regulation)

rownames(plots$regulation) <- name


pw = rownames(selection1)
pw
plots$pw.v.plots$CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS


h2 = Heatmap(all_gsva_tvn_cut[which(rownames(all_gsva_tvn_cut)== "DNA_REPLICATION_INITIATION"),],
        row_names_gp =  grid::gpar(fontsize=7),
        column_names_gp = grid::gpar(fontsize = 7),
        row_names_side = "left",
        show_row_dend = F,
        heatmap_legend_param = list(title = "relative expression"),
        row_names_max_width = unit(14, "cm"),
        show_column_names = T
)

h2
