#install.packages("ggplot2")
library(ggplot2) # wichtig!!
#library(pheatmap)
library(scales)
#install.packages("gridExtra")               # Install gridExtra package
library("gridExtra")
#install.packages("gtools")                  #gtools package used for foldchange & logratio2foldchange function
library(gtools)
#------------------------------


#loading data

#tcga_exp <- readRDS("./data/tcga_tumor_log2TPM.RDS")
#tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")
tumor_vs_norm <- readRDS("./data/tcga_tumor_normal_datascience_proj_2022.RDS")
#genesets <- readRDS("./data/hallmarks_genesets.rds")

#-----------alle plots (am ende des durchlaufs die "#" entfernen--------------------------------------
#med.tn                   # Meian aller exprimierten Gene im normalen Gewebe gegen Median aller expr. Gene im tumor Gewebe
#med.tn.diff              # differenz median.tumor-median.normal

# [...]






#-----------reorganizing and cleaning tumor_vs_norm---------------------

#reorganize (Speicherung aller dataframes der LUAD Patienten in separate dataframes)

luad = tumor_vs_norm[["LUAD"]] # im vektor sind nur die infos zu "LUAD" gespeichert (siehe liste)
luad.tumor = luad[["tumor"]]
luad.normal = luad[["normal"]]
luad.anno = luad[["clinical"]]

luad.tumor = as.data.frame(t(luad.tumor)) # transponieren der dataframe, sodass (cols = Gene) und (rows = Patienten)
luad.normal = as.data.frame(t(luad.normal)) # "-"

#extracting gene IDs and ensemble IDs
tumor_IDs = c()
for (element in colnames(luad.tumor)){
  tumor_IDs = c(tumor_IDs, strsplit(element, "|", fixed = T)[[1]][2])  
}
# -> alle elemente in Colnames(also: Colnames = Gennamen, [class = string]) in luad.tumor werden an der Stelle gespalten, wo sich ein "|" befindet
#... vom geschnittenen string nimmt man [[1]] bedeutet man nimmt das erste "|" und davon [2] also das Zweite Element...
#:  ENSG00000167578.16|RAB4B wird gespalten -> "ENSG00000167578.16" und "RAB4B" davon das 2. Element = RAB4B
#Die erhaltene ID wird dann in einem neuen Vector tumor_IDs gespeichert...

tumor_ensembles = c()
for(element in colnames(luad.tumor)){
  tumor_ensembles = c(tumor_ensembles, strsplit(element, "|", fixed = T)[[1]][1])
}

normal_IDs = c()
for (element in colnames(luad.normal)){
  normal_IDs = c(normal_IDs, strsplit(element, "|", fixed = T)[[1]][2])
}

normal_ensembles = c()
for(element in colnames(luad.normal)){
  normal_ensembles = c(normal_ensembles, strsplit(element, "|", fixed = T)[[1]][1])
}

# setting col names to gene IDs; enumerate duplicates

colnames(luad.tumor) = make.names(tumor_IDs, unique = T)
colnames(luad.normal) = make.names(normal_IDs, unique = T)

rm(element, normal_ensembles, normal_IDs, tumor_ensembles, tumor_IDs, tumor_vs_norm) # remove not needed varriables

# structure: rows = patients, cols = gen ID´s


#checking NA´s
sum(is.na(luad.normal)) # 0 NA´s
sum(is.na(luad.tumor)) # "-"
sum(is.na(luad.anno)) # 148 NA´s
NA_sources <- sort(apply(luad.anno, 2, function(x){sum(is.na(x))}))
which(NA_sources > 0)
NA_sources = NA_sources[31:37]
NA_sources # NA sources not sign.


rm(luad, NA_sources)


#------------------------descriptive analysis------------------------------
###########################################################################


#1. create new data frame with mean, median and variance information of each gen
mmv_df = data.frame(med.normal = apply(luad.normal, 2, median),
                    med.tumor = apply(luad.tumor, 2, median),
                    mean.normal = colMeans(luad.normal),
                    mean.tumor = colMeans(luad.tumor),
                    var.normal = apply(luad.normal, 2, var),
                    var.tumor = apply(luad.tumor, 2, var)
                    )


mmv_df$diff.med = c(mmv_df$med.tumor-mmv_df$med.normal) # neue spalte namens "diff" wird erstellt, mit differenzen der mediane d. genexpressionen von Tumor_gewebe- Normales_gewebe

View(head(mmv_df))# rows = gen ID´s, cols =conditions





g1 = ggplot(mmv_df, aes(y = var.normal, x = mean.normal))+ 
  geom_point()+
  labs(title = "mean-variance plot of the normal tissue", # einfach beschriften des plots
       y = "variance", x = ("mean"))

g2 = ggplot(mmv_df, aes(y = var.tumor, x = mean.tumor))+
  geom_point()+
  labs(title = "mean-variance plot of the tumor tissue", # einfach beschriften des plots
       y = "variance", x = ("mean"))

grid.arrange(g1, g2, nrow = 1)

# checking varriance of the gens -> mean variance plots
p = c(.5, .9)
# removing all low varriance gens, below p% quantile
x = quantile(mmv_df$var.normal, probs = p[1], names = F)
mmv_df_new = mmv_df[mmv_df$var.normal > x,]
dim(mmv_df_new)# 9812 removed gens
View(head(mmv_df_new))
# -> we can continue to work with the shortened data




paste("sum of gens with a varr. = 0 BEFORE cut:", sum(which(mmv_df$var.normal == 0)), ". sum of gens with a varr. = 0 AFTER cut:", sum(which(mmv_df_new$var.normal == 0)))
paste("sum of gens with a varr. = 0 BEFORE cut:", sum(which(mmv_df$var.tumor == 0)), ". sum of gens with a varr. = 0 AFTER cut:", sum(which(mmv_df_new$var.tumor == 0)))
# we still have gens with varr = 0 in the tumor dataset



#2.1 removing all the gens under 50% quantile from luad.normal and luad . tumor:
#later we also will remove the varr = 0 genes coming from tumor dataframe (those gens will be removed in bot normal and tumor df)
genenames = rownames(mmv_df_new)
head(genenames)
luad.normal.new = luad.normal[,genenames]# new df with removed 50% quantile values
dim(luad.normal.new)#  58 9812

#2.2same cut with luad.tumor
luad.tumor.new = luad.tumor[,genenames]

#2.3 removing all gens with varr = 0 (otherwise shapiro test won´t work)

rows = c(which(mmv_df_new$var.tumor != 0))# all gens with varr =! 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
zero.varr = rownames(mmv_df_new[rows,])
zero.varr

luad.normal.cut = luad.normal.new[,zero.varr]# luad. normal without the varr =0 gens from tumor
dim(luad.normal.cut)#  58 9784

luad.tumor.cut = luad.tumor.new[,zero.varr]# luad. tumor without the varr = 0 gens
dim(luad.tumor.cut)#  58 9784

mmv_df_cut = mmv_df_new[zero.varr,]# mmv_df without the varr = 0 gens
dim(mmv_df_cut)#  9784  7

which(mmv_df_cut$var.normal == 0) # there are no varr= gens anymore
which(mmv_df_cut$var.tumor == 0)



g3 = ggplot(mmv_df_cut, aes(y = var.normal, x = mean.normal))+ 
  geom_point()+
  labs(title = "mean-variance plot of the normal tissue", # einfach beschriften des plots
       y = "variance", x = ("mean"))+
  geom_smooth()

g4 = ggplot(mmv_df_cut, aes(y = var.tumor, x = mean.tumor))+
  geom_point()+
  labs(title = "mean-variance plot of the tumor tissue", # einfach beschriften des plots
       y = "variance", x = ("mean"))+
  geom_smooth()

grid.arrange(g3, g4, nrow = 1)

rm(p, x, mmv_df, luad, high.var.genes, luad.normal, luad.tumor, mmv_df, colnames, gennames, genenames, p, x, zero.varr, rows)

#our new data:
View(luad.anno)
View(head(luad.normal.cut))
View(head(luad.tumor.cut))
View(head(mmv_df_cut))





#-------------------plotting median_tumor vs median_normal----------------------


#med.tn = ggplot(mmv_dff_new, aes( x = med.tumor, y = med.normal))+
 # geom_point(col = "black", alpha = 0.5, size = .5)+
  #geom_smooth(method = "lm")+
  #geom_hline(yintercept = c(1,-1), col = "tan4")+
  #geom_vline(xintercept = c(1,-1), col = "tan4")

#med.tn # Median aller exprimierten Gene im normalen Gewebe gegen Median aller expr. Gene im tumor Gewebe







#---------------vergleich genexpression tumor und normalgewebe mithilfe de log2foldchanges-----------

# #1 test auf normalverteilung der gene: shapiro-wilcon test

#shapiro.test(luad.tumor.new$RAB4B)$p.value # p value of one gen in tumor tissue

#p.values.tumor = c(1: ncol(luad.tumor.cut))# shapiro.test

#for (i in 1: ncol(luad.tumor.cut)){
 # p.values.tumor[i] = shapiro.test(luad.tumor.cut[,i])$p.value
#}

#class(p.values.tumor)#we have the p values and want to add em to the mmv_df_cut
#head(p.values.tumor)


#mmv_df_cut$p.value.tumor = p.values.tumor

#View(head(mmv_df_cut))

# do the same for normal-tissue data:
#p.values.normal = c(1: ncol(luad.normal.cut))# warum funktioniert das nicht???????

#for (i in 1: ncol(luad.normal.cut)){
  #p.values.normal[i] = shapiro.test(luad.normal.cut[,i])$p.value
#}

#class(p.values.normal)#we have the p values and want to add em to the mmv_df_cut
#head(p.values.normal)


#mmv_df_cut$p.value.normal = p.values.normal

#View(head(mmv_df_cut))

# lets see how many gens are normally distributed:
#paste("according to the p-value of the shapiro test done with the gens in the TUMOR tissue there are ",length(which(mmv_df_cut$p.value.tumor < 0.05)), "gens NOT normally distributed, in other words: have a p.value < 5%")
#paste("according to the p-value of the shapiro test done with the gens in the NORMAL tissue there are ",length(which(mmv_df_cut$p.value.normal < 0.05)), "gens NOT normally distributed, in other words: have a p.value < 5%")
# to reduce the data we would have to remove all gens that are not normally distributed and have a p.value < 5%???????????????????????
#???????????????????????????????????????????????
#????????????????????????????????????????????????
#???????????????????????????????????????????????
#???????????????????????????????????????????????
#????????????????????????????????????????????????
#???????????????????????????????????????????????
#???????????????????????????????????????????????
#????????????????????????????????????????????????
#???????????????????????????????????????????????
#???????????????????????????????????????????????
#????????????????????????????????????????????????
#???????????????????????????????????????????????


#--------------------get all gens out with p values < 0,5 (not normally distributed)----------------------

#test = mmv_df_cut[which(mmv_df_cut$p.value.tumor > 0.005), ]# create a df only consisting of the genes with p values > 0,5
#dim(test) #5320 genes   9cols
#test.cut = test[which(test$p.value.normal > 0.005),]
#dim(test.cut) #3963 genes   9cols





# #1 test auf normalverteilung der gene: shapiro-wilcon test

shapiro.test(luad.tumor.cut$RAB4B)$p.value # p value of one gen in tumor tissue

for (i in 1: ncol(luad.tumor.cut)){# p values for tumor tissue
 mmv_df_cut$shapiro.tumor[i] = shapiro.test(luad.tumor.cut[,i])$p.value
}
for (i in 1: ncol(luad.normal.cut)){#p values for normal tissue
  mmv_df_cut$shapiro.normal[i] = shapiro.test(luad.normal.cut[,i])$p.value
}

View(head(mmv_df_cut))

#lets view the # normal distr. genes
p = (.5)
qt = quantile(mmv_df_cut$shapiro.tumor, probs = p, names = F)
qn = quantile(mmv_df_cut$shapiro.normal, probs = p, names = F)

a1 = ggplot(mmv_df_cut, aes(shapiro.tumor))+
  geom_density(color="darkblue", fill="lightblue")+
  labs(title = "p-value density shapiro test tumor tissue", x ="p-values", y = "density")+
  geom_vline(xintercept = 0.05, linetype="dashed")+
  geom_vline(xintercept = qt, col = "red")
  
a2 = ggplot(mmv_df_cut, aes(shapiro.normal))+
  geom_density(color="darkblue", fill="lightblue")+
  labs(title = "p-value density shapiro test normal tissue", x ="p-values", y = "density")+
  geom_vline(xintercept = 0.05, linetype="dashed")+
  geom_vline(xintercept = qn, col = "red")


normality.gens = grid.arrange(a1, a2, nrow = 1) # we see that there are a lot of not normally distributed data (more than 50 %)
normality.gens

#------------------------------to see weather the means differ from each other we do the wilcox.test:--------------------------------
p.wilc = c()
for(i in 1: ncol(luad.normal.cut)){
  p.wilc[i] = c(wilcox.test(luad.tumor.cut[,i],luad.normal.cut[,i], paired = TRUE, alternative = "two.sided" )$p.value)
  }



class(p.wilc)
length(p.wilc)
dim(luad.tumor.cut)
head(p.wilc)

mmv_df_cut$p.wilc = p.wilc #-log10 tranformierte p werte, damit sie auf dem plot leichter zu sehen sind
View(head(mmv_df_cut))


mmv_df_cut$diff.mean = c(mmv_df_cut$mean.tumor-mmv_df_cut$mean.normal) # neue spalte namens "diff.mean" wird erstellt, mit differenzen der mittelwerte d. genexpressionen von Tumor_gewebe- Normales_gewebe
View(head(mmv_df_cut))# rows = gen ID´s, cols =conditions



# et voila: a volcano plot
alpha = 0.05 #punkt, ab dem die p Values signifikant werden
#alpha muss noch angepasst werden, da der test über große datenmenge durchgeführt wurde alpha.corr = alpha * 1/#gene
alpha.corr = -log10(alpha/nrow(mmv_df_cut))

for (i in 1:nrow(mmv_df_cut)){ #setting up expression
  ifelse((   (mmv_df_cut$diff.mean[i] > 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ), mmv_df_cut$expr[i] <- "up",
         ifelse((   (mmv_df_cut$diff.mean[i] < 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ),mmv_df_cut$expr[i] <- "down",
                mmv_df_cut$expr[i] <- "not sign."))
}

alpha.corr
# plot for corrected alpha
volcano1 = ggplot(mmv_df_cut, aes(diff.mean, -log10(p.wilc), col = expr))+ #plot
  geom_point(size = .5)+
  labs(title = "volcano plot of high varicance genes -> corrected alpha", y = "-log10(p.value)", x = "log2(Foldchange)")+
  geom_hline(yintercept = alpha.corr,
             linetype = "dashed")+
  scale_color_manual(values=c("blue", "gray", "red"))

volcano1


# plot for ___NOT___ corrected alpha
for (i in 1:nrow(mmv_df_cut)){ #setting up expression
  ifelse((   (mmv_df_cut$diff.mean[i] > 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha)   ), mmv_df_cut$expr[i] <- "up",
         ifelse((   (mmv_df_cut$diff.mean[i] < 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha)   ),mmv_df_cut$expr[i] <- "down",
                mmv_df_cut$expr[i] <- "not sign."))
}

volcano2 = ggplot(mmv_df_cut, aes(diff.mean, -log10(p.wilc), col = expr))+ #plot
  geom_point(size = .5)+
  labs(title = "volcano plot of high varicance genes -> not corrected alpha ", y = "-log10(p.value)", x = "log2(Foldchange)")+
  geom_hline(yintercept = alpha,
             linetype = "dashed")+
  scale_color_manual(values=c("blue", "gray", "red"))

volcano2
grid.arrange(volcano2, volcano1, nrow = 1)


#an alternative to the volcano plot is the manhattan plot:
#we simply plot the log2 fc against the log 2 mean = A
#A = (log2(x) + log2(y))/2 = log2(xy)*1/2}, where x
#and y are respectively the mean of the two groups being compared
# -> in our case: (mean.tumor + mean.normal)/2
A = c((mmv_df_cut$mean.tumor + mmv_df_cut$mean.normal)/2)
head(A)
mmv_df_cut$ma.plot.mean = A

for (i in 1:nrow(mmv_df_cut)){ #setting up expression (for colour)
  ifelse((   (mmv_df_cut$diff.mean[i] > 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ), mmv_df_cut$expr[i] <- "up",
         ifelse((   (mmv_df_cut$diff.mean[i] < 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ),mmv_df_cut$expr[i] <- "down",
                mmv_df_cut$expr[i] <- "not sign."))
}

alpha.corr

ma.plot1 = ggplot(mmv_df_cut, aes(y = diff.mean, x = ma.plot.mean))+ #plot
  geom_point(size = .5, col = " re")+
  labs(title = "ma plot of high varicance genes", y = "log2(t/n) = Foldchange", x = "log2(tn)/2 = Mean")

ma.plot1

#------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# wilc test with only normally distr. data (extracted out of test.cut)
dim(test.cut)
genenames = rownames(test.cut)
head(genenames)

luad.normal.final = luad.normal.cut[,genenames]# create df with final genes in normal
dim(luad.normal.final)
luad.tumor.final = luad.tumor.cut[,genenames]# same with tumor
dim(luad.tumor.final)
mmv_df_final = mmv_df_cut[genenames,]# same with tumor
dim(luad.tumor.final)
View(head(mmv_df_final))







col = c("blue", "red", "gray", "black")
alpha = 12 #punkt, ab dem die p Values signifikant werden
final.volcano.plot = plot(mmv_df_final$diff.med, -log10(mmv_df_final$p.wilc), pch= 20, # volcano plot final
                          col=ifelse((   (mmv_df_final$diff.med > 0) & (-log10(mmv_df_final$p.wilc) >alpha)   ),col[2],
                                     ifelse((   (mmv_df_final$diff.med < 0) & (-log10(mmv_df_final$p.wilc) >alpha)   ),col[1], col[3])),
                          main = "final volcano plot",
                          xlab = "log2(foldchange)",
                          ylab = "-log10(p values)",
                          abline(h = alpha, col = col[4])
                          )
legend("bottomleft",
       legend=c("Down", "Up", "Not Sign."),
       col= col[1:3],
       pch = 20,
       cex=0.8)

#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------






















#------------------------------------------------------------------------------------------------------------------------------------------------

p = c(.01,.99) # Wahrscheinlichkeiten für die quantile, bitte maximal 2 Werte einsetzen, 1. Wert: unter 50%, 2. Wert über 50%

# hinweis: code ab jetzt bis zum ende des abschnittes ausführen


q = quantile(mmv_df$diff.med, probs = p, names = F) # definition des quantils,  names = FALSE: Vektor besteht nicht aus den Quantilnamen (also ohne 0.25% etc Angaben)
#qg = # alle Werte über dem höheren quantil (q_g_rößer)
#qk = # alle Werte unter dem unteren quantil (q_k_leiner)


qd = c(1:length(mmv_df$diff.med)) # quantil_diff.mederenz vektor mit der länge aller gene, dieser wird nach durchführung
# eines Loopes die information enthalten, ob ein gen über dem oberen quantil, unter dem unteren quantil oder zwischen beiden liegt...

for(i in 1:length(mmv_df$diff.med)){ 
  ifelse(mmv_df$diff.med[i] < q[1], # wenn diff.mederenz der mediane < q[1] (also unteres quantil) dann: ...
         qd[i] <- "drunter", #... steht im qd vector an dieser stelle "drunter". wenn nicht ...
         ifelse(mmv_df$diff.med[i] > q[2], #... dann überprüfen ob diff.med > q[2]. wenn wahr, dann...
                qd[i] <- "drüber", # ..."drüber". wenn nein, dann...
                qd[i] <- "zw")) #... "zwischen"
}
rm(i) # i braucht man ja nicht mehr

qd # qd hat die info über jedes gen gespeichert, ob der wert Über dem oberen quantil, unter dem unteren oder dazwischen liegt...
# Visualisiert:
  
med.tn.diff.med = ggplot(mmv_df,
       aes(y = diff.med, x = 1: length(med.normal), colour = qd))+ # aes = aesthetics, colour of points depends on quantiles
  geom_point(size = .5, shape = 20)+ # add point geometry, change size to 0.5
  labs(title = "(differenz median.tumor-median.normal)", # einfach beschriften des plots
                             y = "differenz", x = ("gene"))+
  geom_text(aes(x= 100, y = -10, label = paste(p[1]*100, "%")), # text an position xy, label = % aller werte unterhalb von unterem quantil
            col = "green")+
  geom_text(aes(x= 100, y = 10, label = paste((1-p[2])*100, "%")), # text an position xy, label = % aller werte oberhalb dem oberen quantil
            col = "red")+
  scale_colour_discrete(name = "Expression", # bearbeiten der legende...
                        labels = c("überextrimiert", "wenig differenz", "unterexprimiert"))+
  ylim(-12, 12) # lim y achse
  
  


med.tn.diff # plot differenz median.tumor-median.normal


#-------------------------------Abschnitt ende--------------------------------------------------





rm(luad.anno, luad.normal.cut, luad.normal.new, luad.tumor.cut, luad.tumor.new, mmv_df_cut, mmv_df_new, test.cut, test, genenames, i, p.values.normal, p.values.tumor,p.wilc, rows, zero.varr, final.volcano.plot, cut.volcano.plot)

#-------------------------------------------------
#volcano plot with median

#wilc test
p.wilc = c()
for(i in 1: ncol(luad.normal.cut)){
  p.wilc[i] = c(wilcox.test(luad.tumor.cut[,i],luad.normal.cut[,i], paired = TRUE, alternative = "two.sided" )$p.value)
}

#paste to df
mmv_df_cut$p.wilc = p.wilc #-log10 tranformierte p werte, damit sie auf dem plot leichter zu sehen sind
View(head(mmv_df_cut))

alpha = 0.05
alpha.corr = -log10(alpha/nrow(mmv_df_cut))

for (i in 1:nrow(mmv_df_cut)){ #setting up expression col
  ifelse((   (mmv_df_cut$diff.med[i] > 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ), mmv_df_cut$expr[i] <- "up",
         ifelse((   (mmv_df_cut$diff.med[i] < 0) & (-log10(mmv_df_cut$p.wilc[i]) > alpha.corr)   ),mmv_df_cut$expr[i] <- "down",
                mmv_df_cut$expr[i] <- "not sign."))
}

# plot for medians
volcano3 = ggplot(mmv_df_cut, aes(diff.med, -log10(p.wilc), col = expr))+ #plot
  geom_point(size = .5)+
  labs(title = "volcano plot of high varicance genes -> corrected alpha", y = "-log10(p.value)", x = "log2(Foldchange) with medans")+
  geom_hline(yintercept = alpha.corr,
             linetype = "dashed")+
  scale_color_manual(values=c("blue", "gray", "red"))

volcano3


#PCA:
  






