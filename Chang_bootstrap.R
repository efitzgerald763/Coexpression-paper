library(tidyverse)
library(doParallel)
library(foreach)

load("Microarray_MDD_Normalized_CR_regressed.RData")

prsCor = function(i, gene, datExpr.local) {
  c <- cor.test(datExpr.local[i,], datExpr.local[gene,], use = 'pairwise.complete.obs')
  dfPrs <- data.frame(Gene = rownames(datExpr.local)[i], R = c$estimate, P = c$p.value)
  return(dfPrs)
}

ACC_cont <- datExpr[, datMeta$REGION == "BA9" & datMeta$GROUP == "Control"] 

dfcorACC_cont = data.frame()

## Control
ACC_cont.master <- ACC_cont
sum.res_cont <- NA

num_cores <- 10 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

bootstrap_iteration <- function(b) {
  clmns <- sort(sample(x = seq(1, 28, 1), size = 28, replace = TRUE))
  ACC_cont <- ACC_cont.master[, c(clmns)]
  dfcorACC_cont <- data.frame()
  
  for (i in 1:nrow(ACC_cont)) {
    dfcorACC_cont <- rbind(dfcorACC_cont, prsCor(i, "ENSG00000149485", ACC_cont))
  }
  
  dfcorACC_cont$FDR <- p.adjust(dfcorACC_cont$P, method = "fdr")
  result <- sum(dfcorACC_cont$R > 0.5 & dfcorACC_cont$FDR < 0.05)
  
  name <- paste(b)
  positive <- list(dfcorACC_cont$Gene[dfcorACC_cont$R > 0.5 & dfcorACC_cont$FDR < 0.05])
  negative <- list(dfcorACC_cont$Gene[dfcorACC_cont$R < -0.5 & dfcorACC_cont$FDR < 0.05])
  
  return(list(result = result, positive = positive, negative = negative))
}

results <- foreach(b = 1:10000, .packages = c("dplyr")) %dopar% {
  bootstrap_iteration(b)
}

stopCluster(cl)

sum.res_cont <- sapply(results, function(x) x$result)
my_big_positive_list <- lapply(results, function(x) x$positive)

write.csv(sum.res_cont,"Cont_BA9_size.csv", row.names = FALSE)

saveRDS(my_big_positive_list, file="Control_FADS1_pos.RData")
