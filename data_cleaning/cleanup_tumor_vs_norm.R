library(ggplot2) # wichtig!!
#library(pheatmap)
library(scales)
install.packages("gridExtra")               # Install gridExtra package
library("gridExtra")
install.packages("gtools")                  #gtools package used for foldchange & logratio2foldchange function
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

rm(element, normal_ensembles, normal_IDs, tumor_ensembles, tumor_IDs, tumor_vs_norm, luad, luad.annot) # remove not needed varriables

# structure: rows = patients, cols = gen ID´s

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


mmv_df$diff = c(mmv_df$med.tumor-mmv_df$med.normal) # neue spalte namens "diff" wird erstellt, mit differenzen der mediane d. genexpressionen von Tumor_gewebe- Normales_gewebe

View(head(mmv_df))# rows = gen ID´s, cols =conditions


#-------------------plotting----------------------


med.tn = ggplot(mmv_df, aes( x = med.tumor, y = med.normal))+
  geom_point(col = "black", alpha = 0.5, size = .5)+
  geom_smooth(method = "lm")+
  geom_hline(yintercept = c(1,-1), col = "tan4")+
  geom_vline(xintercept = c(1,-1), col = "tan4")

med.tn # Median aller exprimierten Gene im normalen Gewebe gegen Median aller expr. Gene im tumor Gewebe
#
#
#
#
#
#
#
#
#
#
#---------------vergleich genexpression tumor und normalgewebe mithilfe der differenz-----------

p = c(.01,.99) # Wahrscheinlichkeiten für die quantile, bitte maximal 2 Werte einsetzen, 1. Wert: unter 50%, 2. Wert über 50%

# hinweis: code ab jetzt bis zum ende des abschnittes ausführen


q = quantile(mmv_df$diff, probs = p, names = F) # definition des quantils,  names = FALSE: Vektor besteht nicht aus den Quantilnamen (also ohne 0.25% etc Angaben)
#qg = # alle Werte über dem höheren quantil (q_g_rößer)
#qk = # alle Werte unter dem unteren quantil (q_k_leiner)


qd = c(1:length(mmv_df$diff)) # quantil_differenz vektor mit der länge aller gene, dieser wird nach durchführung
# eines Loopes die information enthalten, ob ein gen über dem oberen quantil, unter dem unteren quantil oder zwischen beiden liegt...

for(i in 1:length(mmv_df$diff)){ 
  ifelse(mmv_df$diff[i] < q[1], # wenn differenz der mediane < q[1] (also unteres quantil) dann: ...
         qd[i] <- "drunter", #... steht im qd vector an dieser stelle "drunter". wenn nicht ...
         ifelse(mmv_df$diff[i] > q[2], #... dann überprüfen ob diff > q[2]. wenn wahr, dann...
                qd[i] <- "drüber", # ..."drüber". wenn nein, dann...
                qd[i] <- "zw")) #... "zwischen"
}
rm(i) # i braucht man ja nicht mehr

qd # qd hat die info über jedes gen gespeichert, ob der wert Über dem oberen quantil, unter dem unteren oder dazwischen liegt...
# Visualisiert:
  
med.tn.diff = ggplot(mmv_df,
       aes(y = diff, x = 1: length(med.normal), colour = qd))+ # aes = aesthetics, colour of points depends on quantiles
  geom_point(size = .5, shape = 20)+ # add point geometry, change size to 0.5
  labs(title = "(Differenz median.tumor-median.normal)", # einfach beschriften des plots
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


























#----noch nicht beachten-------------------------
#---------------2. calculate foldchange-------------------------------------------
(mmv_df$med.tumor[1]/mmv_df$med.normal[1]) # first gene -> 0.9037896 unsignificant
#lets try with MYC
rownumber = which(rownames(mmv_df)== "MYC")
mmv_df$med.tumor[rownumber]/mmv_df$med.normal[rownumber]
# for all rows, add a new col with fc values
fc = (mmv_df$med.tumor/mmv_df$med.normal)
mmv_df$fc = fc
# plot fold changes

ggplot(mmv_df, aes(y = fc, x = 1:19624))+
  geom_point()+
  labs(y = "foldchange", x = "gens")



#----------------------test foldchange werte mit unterschiedlichen methoden-------------------
fct = data.frame(normal = c(2,4,-2,4,-4,-2),
                 tumor = c(4,2,4,-4,-2,-4))
rownames(fct) = make.names(c("Gen1", "Gen2", "Gen3", "Gen4", "Gen5", "Gen6"))
normal.fc = fct$tumor/fct$normal
diff = abs(fct$tumor-fct$normal)/fct$normal


fct$fc = normal.fc
fct$diff.fc = diff
fct$differenz = c(fct$tumor-fct$normal)
fct
#-------------------------------------------------------