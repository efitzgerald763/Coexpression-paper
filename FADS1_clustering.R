library(Seurat)
library(tidyverse)
library(data.table)
library(biomaRt)
library(NatParksPalettes)
library(EnhancedVolcano)
library(ggthemes)
library(ggh4x)
library(speckle)
library(limma)

Astro_sub <- read_rds("Astro_SCT_integ_male_v2.rds")

# Remove batches with too few cells to integrate
Idents(object = Astro_sub) <- "Batch"
Astro_sub <- subset(x = Astro_sub,  idents = c('1M', '2M'), invert = TRUE)
Astro_sub$Condition <- gsub("Case","MDD",Astro_sub$Condition)

DefaultAssay(Astro_sub) <- 'RNA'
Astro_sub[['integrated']]  <- NULL
Astro_sub <- NormalizeData(Astro_sub)

subset_FADS1 <- Astro_sub

#-------------------------------- Subset to FADS1 DLPFC network -------------------------------------
FADS1_DLPFC <- read.table("ENSG00000149485_Brain_Frontal_Cortex_BA9")
FADS1_DLPFC <- FADS1_DLPFC[FADS1_DLPFC$R > 0.5 & FADS1_DLPFC$FDR <0.05,] %>%
  dplyr::select(Gene)
FADS1_DLPFC <- FADS1_DLPFC$Gene
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
FADS1_DLPFC <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                               "external_gene_name", "description"),
                     values=FADS1_DLPFC,mart= mart)
FADS1_DLPFC <- FADS1_DLPFC$external_gene_name
subset_FADS1 <- subset_FADS1[rownames(subset_FADS1) %in% FADS1_DLPFC,]

#----------------------------------------------------------------------------------------------------

subset_FADS1.list <- SplitObject(subset_FADS1, split.by = "Batch")

subset_FADS1.list <- lapply(X = subset_FADS1.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = subset_FADS1.list)

anchors <- FindIntegrationAnchors(object.list = subset_FADS1.list, 
                                  anchor.features = features)

subset_FADS1 <- IntegrateData(anchorset = anchors)

DefaultAssay(subset_FADS1) <- "integrated"

subset_FADS1 <- ScaleData(subset_FADS1, verbose = FALSE)
subset_FADS1 <- RunPCA(subset_FADS1, npcs = 30, verbose = FALSE)
ElbowPlot(object = subset_FADS1, 
          ndims = 30)

subset_FADS1 <- RunUMAP(subset_FADS1, reduction = "pca", dims = 1:10)
subset_FADS1 <- FindNeighbors(subset_FADS1, reduction = "pca", dims = 1:10)
subset_FADS1 <- FindClusters(subset_FADS1, resolution = 0.2)

subset_FADS1$seurat_clusters <- gsub("0","Cluster 0", subset_FADS1$seurat_clusters)
subset_FADS1$seurat_clusters <- gsub("1","Cluster 1", subset_FADS1$seurat_clusters)
subset_FADS1$seurat_clusters <- gsub("2","Cluster 2", subset_FADS1$seurat_clusters)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.7, "cm")
)

DimPlot(subset_FADS1, reduction = "umap", group.by='seurat_clusters') + 
  theme(plot.title = element_blank(),
        legend.position = c(0.77, 0.15),
        axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(x = axis, y = axis) +
  scale_color_manual(values=natparks.pals("Yellowstone"))



# Extract metadata
meta <-subset_FADS1@meta.data

# Rename columns so they're more intuitive
meta <- meta %>%
  dplyr::rename(FADS1_clusters = seurat_clusters)

Astro_sub@meta.data <- meta

# Stacked + percent
ggplot(meta, aes(fill=FADS1_clusters, x=Condition)) + 
  geom_bar(position="fill", stat="count") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=17),
        legend.title = element_blank()) +
  xlab("") +
  ylab("Proportion of astrocytes") +
  scale_fill_manual(values=natparks.pals("Yellowstone"))

#-----------------------------------------------------------------------

Idents(Astro_sub) <- "FADS1_clusters"
markers <- FindAllMarkers(Astro_sub, only.pos = TRUE, min.pct = 0, 
                          logfc.threshold = 0.1)

write.csv(markers, "cluster_markers.csv")
