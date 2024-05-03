library(edgeR)
library(biomaRt)
library(nichenetr)

# Get FADS1 neg genes i.e. geneset of interest
Neg_net <- readRDS("my_big_negative_GTEx_DLPFC_list.RData")
Neg_net <- Neg_net$ENSG00000149485_ %>%
  as.data.frame()
Neg_net <- Neg_net$c..ENSG00000187961....ENSG00000236423....ENSG00000097021....ENSG00000116329...
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
geneset_oi <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                              "external_gene_name", "description"),
                    values=Neg_net,mart= mart)
geneset_oi <- geneset_oi$external_gene_name

# Get FADS1 pos genes i.e. ligands of interest
Pos_net <- readRDS("my_big_positive_GTEx_DLPFC_list.RData")
Pos_net <- Pos_net$ENSG00000149485_ %>%
  as.data.frame()
Pos_net <- Pos_net$c..ENSG00000179403....ENSG00000177133....ENSG00000142611....ENSG00000049245...
FADS_pos <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                           "external_gene_name", "description"),
                 values=Pos_net,mart= mart)
FADS_pos <- FADS_pos$external_gene_name

# Get astrocyte expressed genes
load("Common_astro_clusters.RData")
AvgExp <- AverageExpression(Astro_sub, 
                            group.by = c("Sample"),
                            assays = 'RNA',
                            slot = "counts",
                            return.seurat = FALSE)
AvgExp <- AvgExp$RNA
genes_to_keep <- apply(cpm(AvgExp) > 0.1, 1, sum) > 0.25 * ncol(AvgExp) 
astro_back <- rownames(AvgExp[genes_to_keep, ])

# Get excitatory neuron expressed genes
soExcit <- readRDS("Excitatory_Cont_MDD.rds")
AvgExp <- AverageExpression(soExcit, 
                            group.by = c("Sample"),
                            assays = 'RNA',
                            slot = "counts",
                            return.seurat = FALSE)
AvgExp <- AvgExp$RNA
genes_to_keep <- apply(cpm(AvgExp) > 0.1, 1, sum) > 0.25 * ncol(AvgExp) 
excit_back <- rownames(AvgExp[genes_to_keep, ])


# Get background genes (i.e. those tested for coexp in GTEx)
load("GTEx_DLPFC_regressed.RData")
background <- rownames(datExpr.AllRegressed) 

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
background <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                               "external_gene_name", "description"),
                     values=background,mart= mart)
background <- background$external_gene_name


# Increase timeout threshold
options(timeout=600)

# Load precalulated NicheNet data
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

# Get ligands and receptors in the resource
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

# only keep the intersect between the resource and the data
expressed_ligands <- intersect(ligands, astro_back)
expressed_receptors <- intersect(receptors, excit_back)

# filter the network to only include ligands for which both the ligand and receptor are expressed
potential_ligands <- lr_network %>% 
  filter(from %in% FADS_pos & to %in% expressed_receptors) %>%
  pull(from) %>% unique()


ligand_activities <- predict_ligand_activities(geneset = geneset_oi, 
                                               background_expressed_genes = background,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)



ligand_activities <- ligand_activities %>% 
  arrange(-aupr) %>% 
  mutate(rank = rank(desc(aupr)))

top_ligands <- ligand_activities %>%
  top_n(15, aupr) %>% 
  arrange(-aupr) %>%
  pull(test_ligand) %>%
  unique()


# get regulatory potentials
ligand_target_potential <- map(top_ligands,
                               ~get_weighted_ligand_target_links(.x,
                                                                 geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix,
                                                                 n = 500)) %>%
  bind_rows() %>% 
  drop_na()

# prep for visualization
active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = ligand_target_potential, 
                                      ligand_target_matrix = ligand_target_matrix)

# order ligands & targets
order_ligands <- intersect(top_ligands,
                           colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets <- ligand_target_potential$target %>%
  unique() %>% 
  intersect(rownames(active_ligand_target_links)) %>%
  make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>%
  make.names() 
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>%
  make.names()

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>%
  t()

vis_ligand_target <- vis_ligand_target %>%
  as.data.frame() %>%
  rownames_to_column("ligand") %>%
  as_tibble()

# Convert dot to underscore and set ligand as index
vis_ligand_target$ligand <- gsub("\\.", "_", vis_ligand_target$ligand)
rownames(vis_ligand_target) <- vis_ligand_target$ligand
vis_ligand_target <- vis_ligand_target[, colSums(vis_ligand_target >= 0.05) > 0]

# Plot heatmap
library(ggplot2)
library(reshape2)

melted_data <- melt(vis_ligand_target)

plotme <- ggplot(data = melted_data, aes(x = variable, y = ligand, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "purple", limits = c(0, 0.16)) +
  theme_bw() +
  labs(x = "Prioritized receptors in the FADS1- network", 
       y = "Prioritized ligands in the FADS1+ network", 
       fill = "Regulatory potential") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 16))







