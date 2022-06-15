library(parallel)



#enrichment funktion
enrichment = function(expressiondata, genesets = genesets_ids){
  ESmatrix = sapply(genesets, FUN = function(x){
    ins = na.omit(which(rownames(expressiondata) %in% x == TRUE))#indices der gene im aktuellen set
    outs = -ins#indices der gene nicht im aktuellen set
    #gibt einen vektor der für jeden patienten den pval für das aktuelle gene enthält
    res = NULL
    for (i in 1:ncol(expressiondata)){#testet für jeden patienten
      res[i] = wilcox.test(expressiondata[ins,i],expressiondata[outs,i],'two.sided')$p.value
    }
    return(res)
  })
  row.names(ESmatrix) = colnames(expressiondata); return(ESmatrix)
}
pvalueslist_metabol = mclapply(tumor_type_dfs, FUN = function(x){enrichment(x,metabolism_list_gs)}, mc.cores = 6)


pvaluelist_combined <- list()
for(type in names(pvalueslist)){
  pvaluelist_combined[[type]] <- as.data.frame(cbind(pvalueslist[[type]], pvalueslist_metabol[[type]]))
}




