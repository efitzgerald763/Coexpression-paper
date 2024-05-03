library(tidyverse)
library(data.table)


# Read raw, un-normalized gene read counts from GTEx
gtex.read <- as.data.frame(fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
                                 header = TRUE))

# Read GTEx sample- and subject-level information
gtex.sample.attrib <- as.data.frame(fread("GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", 
                                          header = TRUE))
gtex.pheno <- as.data.frame(fread("GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt", 
                                  header = TRUE))

# Extract attributes for just the brain (based on the "regions" list we previously created)
gtex.brain.attrib <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% 'Brain - Frontal Cortex (BA9)'), ] 


tissue.dat <- read.table("Brain_Frontal_Cortex_BA9.v8.normalized_expression.bed.gz",
                         comment.char = "",header = TRUE)

region.sampids <- gtex.sample.attrib$SAMPID[which(gtex.sample.attrib$SMTSD %in% 'Brain - Frontal Cortex (BA9)')] 
region.samp.att <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% 'Brain - Frontal Cortex (BA9)'), ] 

region.counts <- gtex.read[, which(colnames(gtex.read) %in% region.sampids)]
rownames(region.counts) <- gtex.read$Name 

tpm.use.genes <- gtex.read$Name[which(gtex.read$Name %in% tissue.dat$gene_id)] 

gtex.use.samples <- colnames(tissue.dat)[grep("GTEX", colnames(tissue.dat))] 
gtex.read.refs <- gsub("-", ".", colnames(region.counts)) 
use.samp.idxs <- unlist(lapply(gtex.use.samples, grep, x = gtex.read.refs)) 


datExpr.counts <- region.counts[which(rownames(region.counts) %in% tpm.use.genes), use.samp.idxs]

# Get region-specific covariates
tissue.covars.tmp <- read.table("Brain_Frontal_Cortex_BA9.v8.covariates.txt", comment.char = "", header = TRUE)  
tissue.covar <- tissue.covars.tmp[, grep("GTEX", colnames(tissue.covars.tmp))] 
rownames(tissue.covar) <- tissue.covars.tmp$ID 

tissue.covar <- tissue.covar[(1:18), ] 
region <- rep('Brain_Frontal_Cortex_BA9', dim(tissue.covar)[2]) 

# Phenotype information for the current tissue
region.pheno <- gtex.pheno[gsub("-", ".", gtex.pheno$SUBJID) %in% colnames(tissue.covar), ]

# Sequencing metrics
seq <- region.samp.att[which(region.samp.att$SAMPID %in% colnames(datExpr.counts)), ] 


print(length(which(gsub("-", ".", region.pheno$SUBJID) == colnames(tissue.covar))))
count.sub.names <- unlist(lapply(colnames(datExpr.counts), function(x) paste(strsplit(x, '-')[[1]][1:2], collapse = "."))) 
length(which(count.sub.names == colnames(tissue.covar)))
seq <- seq[match(count.sub.names, colnames(tissue.covar)), ]

t.tissue.covar <- as.data.frame(t(tissue.covar))

sub_num <- rownames(t.tissue.covar)
datMeta <- cbind(t.tissue.covar, sub_num, region, region.pheno, seq) 

  
rownames(datExpr.counts) <- gsub("\\..*", "", rownames(datExpr.counts))



gtex.rsem <- as.data.frame(fread("GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz", 
                                 header = TRUE)) 
gtex.rsem[1:5,1:5]

tx2gene <- data.frame(gtex.rsem[,1:2]) 
rownames(gtex.rsem) <- gtex.rsem$transcript_id 
gtex.rsem <- gtex.rsem[, c(-1,-2)]
gtex.rsem <- as.matrix(gtex.rsem)
gtex.rsem <- rowsum(gtex.rsem, tx2gene$gene_id)
rownames(gtex.rsem) <- gsub("\\..*", "", rownames(gtex.rsem))
gtex.rsem <- gtex.rsem[match(rownames(datExpr.counts), rownames(gtex.rsem)), ]
gtex.rsem <- gtex.rsem[, match(colnames(datExpr.counts), colnames(gtex.rsem))]

save(gtex.rsem, datMeta, file = "GTEX_DLPFC_datExpr_counts_brain_RSEM.RData")
