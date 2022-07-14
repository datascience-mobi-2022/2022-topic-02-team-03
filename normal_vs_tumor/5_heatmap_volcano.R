
library(ComplexHeatmap)
library(pheatmap)

#GSVA c5 results -> volcano & heatmap

all_gsva_tvn <- readRDS("./data/all_gsva_tvn.RDS")
sd = c(apply(all_gsva_tvn, 1, sd))
q = quantile(sd, probs = c(.8,.9,.95,.98,1), names = F)
q

all_gsva_tvn_cut = all_gsva_tvn[which(sd >= q[3]),]

dim(all_gsva_tvn_cut)

Heatmap(all_gsva_tvn_cut,
        row_names_gp =  grid::gpar(fontsize=7),
        row_names_side = "left",
        show_row_dend = F,
        heatmap_legend_param = list(title = "relative expression"),
        row_names_max_width = unit(14, "cm"),
        show_column_names = F
        )
        

# volcano plots
## 1. calculate p values
pw = dim(all_gsva_tvn_cut)[1]

all_gsva_n = all_gsva_tvn[,1:58]
all_gsva_t = all_gsva_tvn[,59:116]

p.wilc = c(1:pw)
for(i in 1:pw){
  p.wilc[i] = c(wilcox.test(all_gsva_t[i,],all_gsva_n[i,], paired = TRUE, alternative = "two.sided", exact = F, correct = F )$p.value)
}

mean_n = c(1:pw)
for(i in 1:pw){
  mean_n[i] = mean(all_gsva_n[i,])
}

mean_t = c(1:pw)
for(i in 1:pw){
  mean_t[i] = mean(all_gsva_t[i,])
}

fc = c(1:pw)
for(i in 1:pw){
  fc[i] = mean_t[i] - mean_n[i]
}

volc_df = data.frame(p.wilc = p.wilc,
                     mean_n = mean_n,
                     mean_t = mean_t,
                     fc = (fc))

rownames(volc_df) = rownames(all_gsva_tvn_cut)

alpha = 0.5
alpha.corr = -log10(alpha/nrow(volc_df))
alpha.corr


# plot for ___NOT___ corrected alpha
for (i in 1:nrow(volc_df)){ #setting up expression
  ifelse((   (volc_df$fc[i] > 0) & (-log10(volc_df$p.wilc[i]) > alpha.corr)   ), volc_df$expr[i] <- "up",
         ifelse((   (volc_df$fc[i] < 0) & (-log10(volc_df$p.wilc[i]) > alpha.corr)   ),volc_df$expr[i] <- "down",
                volc_df$expr[i] <- "not sign."))
}

selection1 = volc_df[which(abs(volc_df$fc) > (.35)),]

volcano = ggplot(volc_df, aes(fc, -log10(p.wilc), col = expr))+ 
  geom_point(size = 4, alpha = .5)+
  labs(title = "pathway expression", y = "-log10(p.value)", x = "log2(Foldchange)")+
  geom_hline(yintercept = alpha.corr,
             linetype = "dashed")+
  scale_color_manual(values=c("blue", "gray48", "red"))+
  geom_text_repel(data = selection1,
                  aes(fc, -log10(p.wilc), label = rownames(selection1)), size = 3.5, hjust = -.01, col = "black")+
  geom_point(data = selection1,
                  aes(fc, -log10(p.wilc)), size = 4, col = "green", alpha = .5)

volcano + theme_light()

