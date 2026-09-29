# creating seurat object for cellchat ----
# following 020 & 030
analysis_step <- "310_seu_for_cellchat"

## load packages
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
library(forcats)
source("Rscripts/helpers.R")

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# create seurat obj for cellchat ----
## Load seurat objects, all cells and annotated cells
seu_all <- readRDS("RDSfiles/seu_010_all.RDS")

epi <- readRDS("RDSfiles/seu_020_epi.RDS")
epi_names <- Cells(epi)

str <- readRDS("RDSfiles/seu_030_str.RDS")
str_names <- Cells(str)

fib <- subset(str, subset = celltype2 == "Fibro")
fib_names <- Cells(fib)
ec <- subset(str, subset = celltype2 == "EC")
ec_names <- Cells(ec)

seu <- subset(seu_all, cells = c(epi_names, fib_names, ec_names))

## Add annotation
### celltype2
seu$celltype2 <- ""
seu$celltype2[epi_names] <- as.character(epi$celltype2)
seu$celltype2[fib_names] <- "Fibro"
seu$celltype2[ec_names] <- "EC"
epi_levels <- levels(epi$celltype2)
seu$celltype2 <- factor(seu$celltype2, levels = c(epi_levels, "Fibro", "EC"))

## subset to pri_T, remove non-tumor epithelial cells and Myocytes
seu <- subset(seu, subset = tissue_type == "T")
seu <- subset(seu, subset = celltype2 %in% c("Non-T", "EEC", "Pariet", "IM"), invert = TRUE)
seu$gsm <- droplevels(seu$gsm)
seu$celltype2 <- droplevels(seu$celltype2)
table(seu$gsm, seu$celltype2)

## remove samples where any of the cell types has < 10 cells
tab <- table(seu$gsm, seu$celltype2)
keep <- rownames(tab)[apply(tab, 1, function(x) all(x >= 10))]
seu <- subset(seu, subset = gsm %in% keep)
seu$gsm <- droplevels(seu$gsm)
table(seu$gsm, seu$celltype2, seu$layer)

# save minimal seurat object ----
seu <- DietSeurat(seu)
saveRDS(seu, "RDSfiles/seu_310_T_epi_ec_fib.RDS")

# clustering ----
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000)
hvg <- VariableFeatures(seu)
# cc.genes <- Seurat::cc.genes
# seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
npcs <- 20
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, dims = 1:npcs)

DimPlot(seu, group.by = "celltype2", cols = "polychrome", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# CellChat analysis per sample ----
load("RDSfiles/CellChatDB.new.RData")
CellChatDB.use <- subsetDB(CellChatDB.new, search = c("Secreted Signaling"), key = "annotation")
seu_all <- readRDS("RDSfiles/seu_310_T_epi_ec_fib.RDS")
seu_all <- NormalizeData(seu_all) # need normalized count for cellchat
gsm <- levels(seu_all$gsm)
cellchat_list <- list()

for (pid in gsm){
  cat("Processing:", pid, "\n")
  seu <- subset(seu_all, subset = gsm == pid)
  cellchat <- createCellChat(object = seu, group.by = "celltype2", assay = "RNA")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat_list[[pid]] <- cellchat
}

names(cellchat_list) <- gsm

comm_list <- lapply(names(cellchat_list), function(pid) {
  df <- subsetCommunication(cellchat_list[[pid]])
  df$gsm <- pid
  df
})

comm_df <- do.call(rbind, comm_list)
openxlsx2::write_xlsx(comm_df, file = file.path(res_path, "comm_df.xlsx"))

# surmmarize per superficial or deep ----
meta <- readRDS("RDSfiles/metadata.RDS")
gsm_s <- meta %>% dplyr::filter(layer == "S") %>% pull(gsm) %>% intersect(gsm)
gsm_d <- meta %>% dplyr::filter(layer == "D") %>% pull(gsm) %>% intersect(gsm)

comm_summary_s <- comm_df %>%
  dplyr::filter(gsm %in% gsm_s) %>% 
  group_by(source, target, ligand, receptor, pathway_name) %>%
  summarise(
    n_patients = n_distinct(gsm),
    mean_prob = mean(prob, na.rm = TRUE),
    median_prob = median(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), desc(mean_prob))
openxlsx2::write_xlsx(comm_summary_s, file = file.path(res_path, "comm_summary_s.xlsx"))

comm_summary_d <- comm_df %>%
  dplyr::filter(gsm %in% gsm_d) %>% 
  group_by(source, target, ligand, receptor, pathway_name) %>%
  summarise(
    n_patients = n_distinct(gsm),
    mean_prob = mean(prob, na.rm = TRUE),
    median_prob = median(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), desc(mean_prob))
openxlsx2::write_xlsx(comm_summary_d, file = file.path(res_path, "comm_summary_d.xlsx"))

saveRDS(cellchat_list, file.path("RDSfiles", "cellchat_list_321.RDS"))

# visualization ----
## deep layer
df <- openxlsx2::read_xlsx("result/310_seu_for_cellchat/comm_summary_d.xlsx")
## Tumor to EC
df_plot <- df %>%
  filter(source == "Tumor", target == "EC") %>%
  group_by(ligand) %>%
  arrange(desc(n_patients), desc(mean_prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lr_label = paste0(ligand, " - ", receptor),
    pathway_name = factor(pathway_name)
  ) %>%
  arrange(desc(n_patients), desc(mean_prob)) %>%
  mutate(
    highlight = ifelse(ligand == "LRG1" & receptor == "ENG",
                       "LRG1-ENG", "Other")
  )
df_plot <- df_plot %>% slice_head(n=10)

ggplot(df_plot, aes(x = fct_reorder(lr_label, n_patients), y = n_patients, fill = highlight)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Other" = "grey35", "LRG1-ENG" = "firebrick")) +
  scale_y_continuous(
    breaks = seq(0, max(df_plot$n_patients) + 1, by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "n patients"
  ) +
  theme_panel() & NoLegend()
ggsave("Tumor_EC.pdf", path = plot_path, width = 40, height = 30, units = "mm")

## Fibro to Tumor
df_plot <- df %>%
  filter(source == "Fibro", target == "Tumor") %>%
  group_by(ligand) %>%
  arrange(desc(n_patients), desc(mean_prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lr_label = paste0(ligand, " - ", receptor),
    pathway_name = factor(pathway_name)
  ) %>%
  arrange(desc(n_patients), desc(mean_prob)) %>%
  mutate(
    highlight = ifelse(ligand == "LRG1" & receptor == "ENG",
                       "LRG1-ENG", "Other")
  )
df_plot <- df_plot %>% slice_head(n=10)

ggplot(df_plot, aes(x = fct_reorder(lr_label, n_patients), y = n_patients, fill = highlight)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Other" = "grey35", "LRG1-ENG" = "firebrick")) +
  scale_y_continuous(
    breaks = seq(0, max(df_plot$n_patients) + 1, by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "n patients"
  ) +
  theme_panel() & NoLegend()
ggsave("Fib_Tum.pdf", path = plot_path, width = 40, height = 30, units = "mm")

## superficial layer >>> not many interactions hit
df <- openxlsx2::read_xlsx("result/310_seu_for_cellchat/comm_summary_s.xlsx")
## Tumor to EC
df_plot <- df %>%
  filter(source == "Tumor", target == "EC") %>%
  group_by(ligand) %>%
  arrange(desc(n_patients), desc(mean_prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lr_label = paste0(ligand, " - ", receptor),
    pathway_name = factor(pathway_name)
  ) %>%
  arrange(desc(n_patients), desc(mean_prob)) %>%
  mutate(
    highlight = ifelse(ligand == "LRG1" & receptor == "ENG",
                       "LRG1-ENG", "Other")
  )
df_plot <- df_plot %>% slice_head(n=10)

ggplot(df_plot, aes(x = fct_reorder(lr_label, n_patients), y = n_patients, fill = highlight)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Other" = "grey35", "LRG1-ENG" = "firebrick")) +
  scale_y_continuous(
    breaks = seq(0, max(df_plot$n_patients) + 1, by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "n patients"
  ) +
  theme_panel() & NoLegend()
ggsave("superficial_Tumor_EC.pdf", path = plot_path, width = 40, height = 30, units = "mm")

## Fibro to Tumor
df_plot <- df %>%
  filter(source == "Fibro", target == "Tumor") %>%
  group_by(ligand) %>%
  arrange(desc(n_patients), desc(mean_prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lr_label = paste0(ligand, " - ", receptor),
    pathway_name = factor(pathway_name)
  ) %>%
  arrange(desc(n_patients), desc(mean_prob)) %>%
  mutate(
    highlight = ifelse(ligand == "LRG1" & receptor == "ENG",
                       "LRG1-ENG", "Other")
  )
df_plot <- df_plot %>% slice_head(n=10)

ggplot(df_plot, aes(x = fct_reorder(lr_label, n_patients), y = n_patients, fill = highlight)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Other" = "grey35", "LRG1-ENG" = "firebrick")) +
  scale_y_continuous(
    breaks = seq(0, max(df_plot$n_patients) + 1, by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "n patients"
  ) +
  theme_panel() & NoLegend()
ggsave("superficial_Fib_Tum.pdf", path = plot_path, width = 40, height = 30, units = "mm")