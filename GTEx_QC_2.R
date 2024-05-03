library(tidyverse)
library(edgeR)
library(WGCNA)
library(corrplot)
library(lme4)


load("GTEX_DLPFC_datExpr_counts_brain_RSEM.RData")
datExpr.counts <- gtex.rsem
rm(gtex.rsem)

datMeta <- datMeta %>% 
  filter(!is.na(datMeta$SMRDLGTH), # Maximum read length
         INCEXC == "TRUE", # 	Eligible For Study
         MHALS != "1", # Amyotropic Lateral Sclerosis
         MHALZDMT != "1", #Alzheimer's OR Dementia
         MHDMNTIA != "1", #Dementia With Unknown Cause
         MHENCEPHA != "1", #	Active Encephalitis
         MHFLU != "1", #Influenza (acute viral infection including avian influenza)
         MHMS != "1", #Multiple Sclerosis
         MHPRKNSN != "1", #Parkinson's Disease
         MHREYES != "1", #Reyes Syndrome
         MHSEPSIS != "1", #Documented Sepsis
         MHLUPUS != "1", #Systemic Lupus
         MHCVD != "1", #Cerebrovascular Disease (stroke, TIA, embolism, aneurysm, other circulatory disorder affecting the brain)
         MHALZHMR != "1", #Alzheimer's
         MHBCTINF != "1", #Bacterial Infections (including septicemia (bacteria in the blood), meningococcal disease, staphylococcal infection, streptococcus, sepsis)
         MHCANCERC != "1", #Current Diagnosis Of Cancer
         MHDPRSSN != "1", #MDD
         MHHEROIN != "1", #Heroin Use
         MHPLLABS != "1", #Prescripion pill abuse
         MHPNMIAB != "1", #Pneumonia
         MHSCHZ != "1", #Schizophrenia
         MHSZRSU != "1", #Unexplained seizures
         MHPSBLDCLT != "1") #Positive blood cultures


datExpr.counts <- datExpr.counts[, match(datMeta$SAMPID, colnames(datExpr.counts))] #Gets the gene counts for remaining samples

# Interval of onset to death- Protected data
datMeta$DTHCODD <- as.numeric(datMeta$DTHCODD)
datMeta$DTHCODD[datMeta$DTHCODD >= 0 & datMeta$DTHCODD < 2] <- "0to2h"
datMeta$DTHCODD[datMeta$DTHCODD >= 2 & datMeta$DTHCODD < 10] <- "2to10h"
datMeta$DTHCODD[datMeta$DTHCODD >= 10 & datMeta$DTHCODD < 72] <- "10hto3d"
datMeta$DTHCODD[datMeta$DTHCODD >= 72 & datMeta$DTHCODD < 504] <- "3dto3w"
datMeta$DTHCODD[datMeta$DTHCODD >= 504] <- "3wplus"
datMeta$DTHCODD[is.na(datMeta$DTHCODD)] <- "unknown"

# Filter genes using edgeR package
genes_to_keep <- apply(cpm(datExpr.counts) > 0.1, 1, sum) > 0.25 * ncol(datExpr.counts)
table(genes_to_keep)
# FALSE  TRUE 
# 1350   23326 

# Nomalize
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr.counts[genes_to_keep, ]), method = "TMM"), log = TRUE)
par(mar=c(1,1,1,1))
boxplot(datExpr.norm, range = 0, main = "Filtered, normalized counts", axes=FALSE)

# Remove outliers by Z-score ----
normadj <- (0.5+0.5*bicor(datExpr.norm))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,pch=19, main="Outlier Detection", xlab="", ylab="Standardized Network\nConnectivity (Z score)")
abline(h=-3, lty=2)
excludesampleID <- as.data.frame(names(which(z.ku < -3)))

# 4 outliers to remove
datMeta <- datMeta[-match(t(excludesampleID), datMeta$SAMPID), ] # 107 samples remaining

# Construct PCs from sequencing metrics

seqMet <- c('SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMMAPRT', 'SMEXNCRT',
            'SMGNSDTC', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 
            'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 
            'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 
            'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 
            'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 
            'SMRRNART', 'SME1MPRT', 'SME2PCTS') # 'SMRDLGTH' skipped due to zero variance
datSeq <- datMeta[, seqMet] # 32 columns

seqPCs <- prcomp(t(scale(datSeq, scale=T)), center=T)
topPC <- seqPCs$rotation[,1:10];
varexp <- (seqPCs$sdev)^2 / sum(seqPCs$sdev^2)
topvar <- varexp[1:10]
for(i in 1:10) datMeta[,paste0("seqPC",i)] = topPC[,i]


plot(datMeta$SMRIN, datMeta$SMMPUNRT) # RIN and Mapped Unique Rate of Total: Ratio of mapping of reads that were aligned and were not duplicates to total reads
summary(lm(datMeta$SMMPUNRT ~ datMeta$SMRIN)) # No effect p=0.289


# Reconstruct model matrix after removing outliers
mod1 <- model.matrix(~ seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + 
                             SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + AGE, datMeta)

colnames(mod1) <- make.names(colnames(mod1))

datExpr.norm <- datExpr.norm[, match(datMeta$SAMPID, colnames(datExpr.norm))]

Y <- as.matrix(datExpr.norm)
X <- mod1

# Calculate betas with ordinary least squares 
beta <- (solve(t(X) %*% X) %*% t(X)) %*% t(Y)

# Subtract betas from expression matrix (except intercept)
datExpr.AllRegressed <- Y - t(X %*% beta)

datMeta <- as.data.frame(datMeta)
datMeta$X <- rownames(datMeta)

save(datMeta, datExpr.AllRegressed, file = "GTEx_DLPFC_regressed.RData")

