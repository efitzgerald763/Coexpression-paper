library(tidyverse)

load("GTEx_DLPFC_regressed.RData")

prsCor = function(i, gene, datExpr) {
	c <- cor.test(datExpr[i,], datExpr[gene,], use = 'pairwise.complete.obs')
	dfPrs <- data.frame(Gene = rownames(datExpr)[i], R = c$estimate, P = c$p.value)
	return(dfPrs)
}

data <- read.table("HMAG_genes.txt", header=FALSE) %>%
  dplyr::select(V2)

my_big_positive_list <- list()
my_big_negative_list <- list()
All_networks <- list()

for(r in 1:nrow(data)) {
  print(data[r,])
  try({ 
    dfcor <- data.frame()
    for (i in 1:nrow(datExpr.AllRegressed)) {
      if (i%%2000 == 0) {print(i)}
      dfcor <- rbind(dfcor, prsCor(i, data[r,1], datExpr.AllRegressed))
    }
    dfcor$FDR <- p.adjust(dfcor$P, method = "fdr")
    
    positive <- list(dfcor$Gene [dfcor$R > 0.5 & dfcor$FDR < 0.05])
    negative <- list(dfcor$Gene [dfcor$R < -0.5 & dfcor$FDR < 0.05])
    name <- paste(data[r, 1], data[r,2], sep = "_")
    my_big_positive_list[[name]] <- positive
    my_big_negative_list[[name]] <- negative
    All_networks[[name]] <- dfcor
  }, silent = FALSE)	
}

saveRDS(my_big_positive_list, file="my_big_positive_GTEx_DLPFC_list.RData")
saveRDS(my_big_negative_list, file="my_big_negative_GTEx_DLPFC_list.RData")
saveRDS(All_networks, file="All_GTEx_DLPFC_networks.RData")
