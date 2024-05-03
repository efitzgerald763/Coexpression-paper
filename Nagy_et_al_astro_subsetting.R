library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(glmGamPoi)
library(biomaRt)
DLPFC <- Read10X(data.dir = "/Barcodes/",  gene.column = 1)
so <- CreateSeuratObject(counts = DLPFC, project = "DLPFC")

#Add QC info to the metadata
so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
so$mitoRatio <- PercentageFeatureSet(object = so, pattern = "^MT-")

# Extract metadata, ID group and batch from cell names and rename some columns
metadata <- so@meta.data
metadata$cells <- rownames(metadata)

# ID astrocytes
metadata$Ast <- ifelse(grepl("Ast", rownames(metadata)), "Astrocytes", "Other")

# Extract sample from cell names
metadata$Sample <- word(metadata$cells,1,sep = "\\.")

# Combine with sample info from Corina
MET <- openxlsx::read.xlsx("Metadata.xlsx")%>%
  dplyr::select(Sample, Condition, Bank, Batch, Chemistry, Sequencing, Age, Race, Sex)

x <- merge(metadata, MET, by = 'Sample')
rownames(x) <- x$cells

# Rename columns so they're more intuitive
x <- x %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
so@meta.data <- x


# Extract counts
counts <- GetAssayData(object = so, slot = "counts")

# Output a logical vector for every gene with more than 0 counts per cell
nonzero <- counts > 0

# Only keeping those genes expressed in more than 10 cells
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
so <- CreateSeuratObject(filtered_counts, meta.data = so@meta.data)

# Filter out low quality reads using selected thresholds 
so <- subset(so, subset = nGene > 700 & nGene < 7000)


#-----------------------------------------------------------
#-----------------------------------------------------------


# Run SCTransform on each batch separately
split_seurat <- SplitObject(so, split.by = "Batch")
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], method = "glmGamPoi",  verbose = T)
}

features <- SelectIntegrationFeatures(object.list = split_seurat)
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = split_seurat, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)
so.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

so.integrated <- ScaleData(so.integrated, verbose = FALSE)
so.integrated <- RunPCA(so.integrated, verbose = FALSE)
so.integrated <- RunUMAP(so.integrated, dims = 1:50)

DimPlot(so.integrated, group.by = "Batch")
FeaturePlot(so.integrated, features = c("SLC1A2", "GFAP", "GLUL", "AQP4", "SLC1A3"))

Astro_sub <- subset(so.integrated, subset = Ast == 'Astrocytes')

saveRDS(so.integrated, file="Full_SCT_male_v2.rds")
saveRDS(Astro_sub, file="Astro_SCT_integ_male_v2.rds")

