library(paletteer)
library(ggsignif)

#------------- Negative network in DLPFC----------------
so <- readRDS("Astro_obj.rds")



# FADS1 DLPFC pos
DLPFC_FADS1_pos <- readRDS("my_big_positive_GTEx_DLPFC_list.RData")
DLPFC_FADS1_pos <- DLPFC_FADS1_pos$ENSG00000149485_
DLPFC_FADS1_pos <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                               "external_gene_name", "description"),
                     values=DLPFC_FADS1_pos,mart= mart)
DLPFC_FADS1_pos <- DLPFC_FADS1_pos$external_gene_name %>%
  as.data.frame()



# Create score
DefaultAssay(Astrocytes) <- 'RNA'
Astrocytes <- NormalizeData(Astrocytes)
Astrocytes <- AddModuleScore(object = Astrocytes, features = DLPFC_FADS1_pos, name = "FADS1_DLPFC_pos")

Ast_met <- Astrocytes@meta.data

# Mean mod score by sample and condition
mean_fads1 <- aggregate(FADS1_DLPFC_pos1 ~ Sample + Condition, data = Ast_met, FUN = mean)

mean_fads1$Condition <- gsub("Case", "MDD", mean_fads1$Condition)

modscore <- ggplot(mean_fads1, aes(x=Condition, y=FADS1_DLPFC_pos1))+ 
  geom_boxplot(aes(fill=Condition),color="#535352", show.legend = FALSE) +
  geom_point() +
  scale_fill_manual(values=c('pink3', 'pink4')) + 
  geom_signif(comparisons = list(c("Control", "MDD")), 
              map_signif_level = F) +
  theme_bw() +
  xlab('') +
  ylab('Module score') +
  theme(legend.position = 'bottom',
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16))


