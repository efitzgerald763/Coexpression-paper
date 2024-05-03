library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

load("Astro_clusters.RData")

counts <- Astro_sub@assays$RNA
subset_FADS1@assays$RNA <- counts

subset_fad <- subset_FADS1
DefaultAssay(subset_fad) <- 'RNA'

subset_fad <- as.SingleCellExperiment(subset_fad)

reducedDim(subset_fad, "PCA", withDimnames=TRUE) <- subset_FADS1[['pca']]@cell.embeddings
reducedDim(subset_fad, "UMAP", withDimnames=TRUE) <- subset_FADS1[['umap']]@cell.embeddings


subset_fad <- miloR::Milo(subset_fad)


subset_fad <- buildGraph(subset_fad, k = 20, d = 20)
subset_fad <- makeNhoods(subset_fad, prop = 0.2, k = 20, d=20, refined = TRUE)

plotNhoodSizeHist(subset_fad)

subset_fad <- countCells(subset_fad, meta.data = data.frame(colData(subset_fad)), sample="Sample")


astro_design <- data.frame(colData(subset_fad))[,c("Sample", "Age", "Batch", "Condition")]

astro_design$Batch <- as.factor(astro_design$Batch) 
astro_design <- distinct(astro_design)
rownames(astro_design) <- astro_design$Sample

subset_fad <- calcNhoodDistance(subset_fad, d=30)

da_results <- testNhoods(subset_fad, design = ~ Batch + Age + Condition, design.df = astro_design)

subset_fad <- buildNhoodGraph(subset_fad)

plotNhoodGraphDA(subset_fad, da_results, layout="UMAP",alpha=0.1) 


#-----------------------------------------------
da_results <- annotateNhoods(subset_fad, da_results, coldata_col = "seurat_clusters")
head(da_results)

Bswarm_nagy <- plotDAbeeswarm(da_results, group.by = "seurat_clusters") +
  xlab("") +ylab("Log fold change") +
  theme(panel.grid = element_blank())



#----------------------------------------------------
#---------------------- DEGs ------------------------
#-----------------------------------------------------
subset_fad <- logNormCounts(subset_fad)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)

keep.rows <- rowSums(logcounts(subset_fad)) != 0
subset_fad <- subset_fad[keep.rows, ]

da_nhood_markers <- findNhoodGroupMarkers(subset_fad, da_results, subset.row = rownames(subset_fad))

options(ggrepel.max.overlaps = 'Inf')

top <- da_nhood_markers[da_nhood_markers$logFC_1 > 0,]
top <- slice_min(top, n=10, adj.P.Val_1)

bottom <- da_nhood_markers[da_nhood_markers$logFC_1 < 0,]
bottom <- slice_min(bottom, n=10, adj.P.Val_1)

label_me <- c(top$GeneID, bottom$GeneID)

library(EnhancedVolcano)

keyvals <- ifelse(
  da_nhood_markers$logFC_1 < 0 & da_nhood_markers$adj.P.Val_1 < 0.05, 'firebrick',
  ifelse(da_nhood_markers$logFC_1 > 0 & da_nhood_markers$adj.P.Val_1 < 0.05, 'forestgreen',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'firebrick'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'forestgreen'] <- 'low'

Volcano <- EnhancedVolcano(da_nhood_markers,
                           lab = da_nhood_markers$GeneID,
                           x = 'logFC_1',
                           y = 'adj.P.Val_1',
                           colAlpha = 1,
                           title = '',
                           subtitle = "",
                           labSize = 4.0,
                           selectLab = label_me,
                           drawConnectors = TRUE,
                           arrowheads = FALSE,
                           FCcutoff = 0,
                           colCustom = keyvals,
                           ylab = bquote(~-log[10]~ '(Adj P-Value)'),
                           xlab = bquote(~log[2]~ 'fold change'),
                           cutoffLineType = 'blank',
                           cutoffLineCol = 'black',
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE) +
  scale_x_continuous(limits = c(-2,2))  +
  scale_y_continuous(limits = c(0,260))  +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  annotate("text", x = -1, y = 250, size = 14/.pt, label = "Enriched in \n other nhoods",
           colour = "firebrick")+
  annotate("text", x = 1, y = 250, size = 14/.pt,label = "Enriched in \n MDD depleted nhoods",
           colour = "forestgreen")


sig <- da_nhood_markers[da_nhood_markers$adj.P.Val_1 < 0.05,]


Upreg <- sig[sig$logFC_1 > 0.25,]%>%
  dplyr::select(GeneID)

background <- da_nhood_markers$GeneID

Upreg_GO <- enrichGO(gene          = Upreg$GeneID,
                     OrgDb         = "org.Hs.eg.db",
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 1,
                     universe = background)

Upreg_GO_res <- Upreg_GO@result

Upreg_GO_res <- Upreg_GO_res[with(Upreg_GO_res,order(p.adjust)),]
Upreg_GO_res <- Upreg_GO_res[1:10,]
Upreg_GO_res$Cluster <- 'Enriched in MDD depleted nhoods'


Downreg <- sig[sig$logFC_1 < -0.25,]%>%
  dplyr::select(GeneID)
Downreg_GO <- enrichGO(gene          = Downreg$GeneID,
                     OrgDb         = "org.Hs.eg.db",
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 1,
                     universe = background)

Downreg_GO_res <- Downreg_GO@result

Downreg_GO_res <- Downreg_GO_res[with(Downreg_GO_res,order(p.adjust)),]
Downreg_GO_res <- Downreg_GO_res[1:10,]
Downreg_GO_res$Cluster <- 'Enriched in other nhoods'


# Capatilize the first letter of each GO term
Capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

Total <- rbind(Upreg_GO_res, Downreg_GO_res)
Total$Description <- Capitalize(Total$Description)
Total$logPvalue <- -log10(Total$p.adjust)


reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

DEG_GO <- ggplot(Total, aes(reorder_within(Description, logPvalue, Cluster), logPvalue, fill = Cluster)) +
  geom_segment(aes(xend = reorder_within(Description, logPvalue, Cluster), yend = 0) 
  ) +
  geom_point(size = 5, aes(colour = Cluster)) + 
  coord_flip() +
  theme_bw() +
  scale_x_reordered() +
  facet_wrap(facets = "Cluster", scales = "free") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) +
  scale_color_manual(values=c("forestgreen","firebrick"))+ 
  geom_hline(yintercept=1.3, linetype="dashed", color = "grey44") +
  xlab("") +
  ylab(expression(paste("-log"[10]," (adj P-value)"))) 


