library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(viridis)

data_dir <-"Br2720_all/Br2720_ant2/"

brain <-Load10X_Spatial(data.dir = data_dir, filename= 'filtered_feature_bc_matrix.h5')

# Change the co-ordinates from characters to integers
brain@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(brain@images[["slice1"]]@coordinates[["tissue"]])
brain@images[["slice1"]]@coordinates[["row"]] <- as.integer(brain@images[["slice1"]]@coordinates[["row"]])
brain@images[["slice1"]]@coordinates[["col"]] <- as.integer(brain@images[["slice1"]]@coordinates[["col"]])
brain@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(brain@images[["slice1"]]@coordinates[["imagerow"]])
brain@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(brain@images[["slice1"]]@coordinates[["imagecol"]])

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial")

brain <- RunPCA(brain, assay = "SCT")
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# Import cluster markers from Nagy analysis and create module scores
markers <- read.csv("cluster_markers.csv")
cluster0 <- markers[markers$cluster == 'Cluster 0' & markers$p_val_adj < 0.05,]
features.to.score <- cluster0$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Cluster0")
P0 <- SpatialFeaturePlot(brain, features = 'Cluster01', pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")  +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))+
  ggtitle('Cluster 0 Score')

cluster1 <- markers[markers$cluster == 'Cluster 1' & markers$p_val_adj < 0.05,]
features.to.score <- cluster1$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Cluster1")
P1 <- SpatialFeaturePlot(brain, features = 'Cluster11', pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")  +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))+
  ggtitle('Cluster 1 Score')

cluster2 <- markers[markers$cluster == 'Cluster 2' & markers$p_val_adj < 0.05,]
features.to.score <- cluster2$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Cluster2")
P2 <- SpatialFeaturePlot(brain, features = 'Cluster21', pt.size.factor = 3) +
  scale_fill_viridis(option="inferno") +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))+
  ggtitle('Cluster 2 Score')


cowplot::plot_grid(P0,P1,P2, nrow = 1)


#----------- Correlate cell types -----------------

full <- read.csv("Full_Nagy_cluster_markers.csv")

Endo <- full[full$cluster == 'Endo' & full$p_val_adj < 0.05,]
features.to.score <- Endo$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Endo")
SpatialFeaturePlot(brain, features = "Endo1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

OPC <- full[full$cluster == 'OPC' & full$p_val_adj < 0.05,]
features.to.score <- OPC$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "OPC")
SpatialFeaturePlot(brain, features = "OPC1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

Excitatory <- full[full$cluster == 'Excitatory' & full$p_val_adj < 0.05,]
features.to.score <- Excitatory$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Excitatory")
SpatialFeaturePlot(brain, features = "Excitatory1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

Inhib <- full[full$cluster == 'Inhib' & full$p_val_adj < 0.05,]
features.to.score <- Inhib$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Inhib")
SpatialFeaturePlot(brain, features = "Inhib1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

Oligo <- full[full$cluster == 'Oligo' & full$p_val_adj < 0.05,]
features.to.score <- Oligo$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Oligo")
SpatialFeaturePlot(brain, features = "Oligo1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

Micro <- full[full$cluster == 'Micro' & full$p_val_adj < 0.05,]
features.to.score <- Micro$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Micro")
SpatialFeaturePlot(brain, features = "Micro1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")

Astro <- full[full$cluster == 'Astro' & full$p_val_adj < 0.05,]
features.to.score <- Astro$gene %>%
  as.data.frame()
brain <- AddModuleScore(object = brain, features = features.to.score, name = "Astro")
SpatialFeaturePlot(brain, features = "Astro1", pt.size.factor = 3) +
  scale_fill_viridis(option="inferno")


###### Make Corrplot ######
library(corrplot)
meta <- brain@meta.data
meta <- meta%>%
  dplyr::rename('Cluster 0 score'=Cluster01)%>%
  dplyr::rename('Cluster 1 score'=Cluster11)%>%
  dplyr::rename('Cluster 2 score'=Cluster21)%>%
  dplyr::rename('Endothelial score'=Endo1)%>%
  dplyr::rename('Microglial score'=Micro1)%>%
  dplyr::rename('Excitatory neuron score'=Excitatory1)%>%
  dplyr::rename('Inhibitory neuron score'=Inhib1)%>%
  dplyr::rename('Astrocyte score'=Astro1)%>%
  dplyr::rename('Oligodendrocyte score'=Oligo1)%>%
  dplyr::rename('OPC score'=OPC1)

meta <- meta%>%
  dplyr::select('Cluster 0 score', 'Cluster 1 score', 'Cluster 2 score','Endothelial score',
                'Microglial score', 'Excitatory neuron score', 'Inhibitory neuron score',
                'Oligodendrocyte score', 'OPC score')

# Make correlation matrix
correlation_matrix <- round(cor(meta),2)


# Generate correlations for all slices individually, before getting the mean correlation coeffecient across all samples for plotting
