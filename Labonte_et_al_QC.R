library(tidyverse)
library(SummarizedExperiment)
library(edgeR)
library(nlme)
library("XML")
library(GEOquery)

datExpr <- fread("Labonte_GSE102556_HumanMDD_fpkmtab.txt.gz")%>%
  as.data.frame()
rownames(datExpr) <- datExpr$gene_id
colnames(datExpr) <- paste("X", colnames(datExpr), sep = "_")

# Read in metadata and create unique IDs that match the expression matrix
datMeta <- read_csv("Labonte_meta.csv")
datMeta$sampIDs <- paste(datMeta$indiv,datMeta$region,sep=".")
datMeta$sampIDs <- paste("X", datMeta$sampIDs, sep = "_")

# Filter for ROI
datMeta <- datMeta[datMeta$region == 'BA25',]
datExpr <- datExpr[, match(datMeta$sampIDs, colnames(datExpr))]

table(datMeta$region)
datMeta$RIN.squared <- datMeta$RIN^2
datMeta$age.squared <- datMeta$Age^2

datMeta$Gender <- factor(datMeta$Gender, levels = c("Male", "Female"))



# -----------------------------------------------------------

mod1 <- model.matrix(~ Condition + Age + age.squared + pH + Gender + PMI + RIN + RIN.squared, data = datMeta)
colnames(mod1) <- make.names(colnames(mod1))

# Filter genes (more than 0.5 FPKM in greater than 25% of samples)
genes_to_keep <- apply(datExpr > 0.5, 1, sum) > 0.25 * ncol(datExpr)
table(genes_to_keep)
 

datExpr <- datExpr[genes_to_keep, ]

# -----------------------------------------------------------

Y <- as.matrix(datExpr)
X <- mod1
beta <- (solve(t(X) %*% X) %*% t(X)) %*% t(Y)

datExpr.AllRegressed <- Y - t(X[, c(3:ncol(X))] %*% beta[c(3:ncol(X)), ])

save(datMeta, datExpr.AllRegressed, file = "Labonte_BA25_regressed.RData")

