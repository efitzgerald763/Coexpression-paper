library(tidyverse)
library(SummarizedExperiment)
library(edgeR)
library(nlme)

##-----Load expression and meta data
datMeta1 <- read.csv("CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv")
datMeta2 <- read.csv("CMC_MSSM-Penn-Pitt_Clinical.csv")
datMeta3 <- read.delim("CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_GemToolsAncestry.tsv")
datMeta <- merge(datMeta1, datMeta2, by="Individual_ID")
datMeta <- merge(datMeta, datMeta3, by="Genotyping_Sample_ID")

rownames(datMeta) <- datMeta$DLPFC_RNA_Sequencing_Sample_ID.x
datMeta$Batch <- as.factor(datMeta$DLPFC_RNA_Sequencing_Library_Batch)
datMeta$RIN2 <- datMeta$DLPFC_RNA_isolation_RIN^2

Age <- datMeta$Age_of_Death
Age[Age=="90+"] <- 90
Age <- as.numeric(Age)
datMeta$Age <- Age
datMeta$Age.squared <- datMeta$Age^2

datMeta$sex <- factor(datMeta$Gender, levels = c("Male", "Female"))

datMeta <- subset(datMeta, grepl("Control", Dx))

columns_of_interest <- c("Age", "Age.squared", "Batch", "pH", "Gender", "PMI_hrs", "DLPFC_RNA_isolation_RIN", "RIN2")
datMeta <- datMeta[complete.cases(datMeta[, columns_of_interest]), ]

datExpr <- read.table("CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv", skip =1, row.names = 1)
id <- read.table("CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv", nrow=1)
colnames(datExpr) <- id

genes_to_keep <- apply(cpm(datExpr) > 0.1, 1, sum) > 0.25 * ncol(datExpr)
table(genes_to_keep)
# FALSE  TRUE 
# 30033 26599 

datExpr <- as.matrix(datExpr)
datExpr <- datExpr[, match(rownames(datMeta), colnames(datExpr))]

datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr[genes_to_keep, ]), method = "TMM"), log = TRUE)

normadj <- (0.5+0.5*bicor(datExpr.norm))^2 
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,pch=19, main="Outlier Detection", xlab="", ylab="Standardized Network\nConnectivity (Z score)")
abline(h=-3, lty=2)
outlier_sample_ids <- names(which(z.ku < -3))


datMeta <- datMeta[!rownames(datMeta) %in% outlier_sample_ids, ] 
datExpr <- datExpr[, match(rownames(datMeta), colnames(datExpr))]

datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr[genes_to_keep, match(rownames(datMeta),colnames(datExpr))]), 
                                    method = "TMM"), log = TRUE)

datSeq <- datMeta[,c(20:29)]; rownames(datSeq) = datMeta$DLPFC_RNA_Sequencing_Sample_ID.x
picard <- read.csv("PicardToolsQC.csv",row.names = 1)
idx <- match(rownames(datSeq), rownames(picard))
datSeq <- cbind(datSeq, picard[idx,])
a <- apply(is.na(datSeq),1,any) | (rownames(datSeq)=="PITT_RNA_BP_PFC_686") #Duplicate sample
datSeq <- datSeq[!a,]; datMeta=datMeta[!a,]
idx <- which(datSeq[1,]>100)
datSeq[,idx] = log10(datSeq[,idx])
PC <- prcomp(t(scale(datSeq, scale=T)), center=T)
topPC <- PC$rotation[,1:10];
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topvar= varexp[1:10]
for(i in 1:10) datMeta[,paste0("seqPC",i)] = topPC[,i]

datExpr.norm <- datExpr.norm[, match(rownames(datMeta), colnames(datExpr.norm))]

# Construct model matrix after removing outliers
mod1 <- model.matrix(~ Age + Age.squared + Batch + pH + Gender + PMI_hrs + DLPFC_RNA_isolation_RIN + RIN2 + 
                       seqPC1 + seqPC2 + seqPC3 + seqPC4+ seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10, 
                     data = datMeta)
colnames(mod1) <- make.names(colnames(mod1))



Y <- as.matrix(datExpr.norm)
X <- mod1

beta <- (solve(t(X) %*% X) %*% t(X)) %*% t(Y)

datExpr.AllRegressed <- Y - t(X[, c(2:ncol(X))] %*% beta[c(2:ncol(X)), ])

datMeta <- as.data.frame(datMeta)
datMeta$X <- rownames(datMeta)

save(datMeta, datExpr.AllRegressed, file = "CMC_control_regressed.RData")


