# 2026-03-19
# Setting up ----

## Make directories
analysis_step <- "020_epi"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(openxlsx2)
library(Seurat)
library(msigdbr)
library(UCell)
library(ggpubr)
library(pals)
source("Rscripts/helpers.R")

## Helper functions
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 5) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering w/o integration ----
seu <- readRDS("RDSfiles/seu_010_all.RDS")
seu <- subset(seu, subset = celltype1 == "Epi")

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

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "layer", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Layer")
ggsave("layer.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# Check markers ----
# features <- readLines("aux_data/gene_set/global_markers.txt")
# sapply(features, save_fp, seu, fp_path)

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

# re-clustering for visualization
seu <- subset(seu, subset = seurat_clusters %in% c(13,15,16), invert = TRUE) # remove obvious contamination
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000)
hvg <- VariableFeatures(seu)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
npcs <- 20
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, dims = 1:npcs)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

fp_path2 <- paste0(fp_path, "2")
fs::dir_create(fp_path2)
add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path2)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers2.csv"))

## annotation
seu$celltype2 <- "Tumor"
seu$celltype2[seu$seurat_clusters %in% c(3,5,10)] <- "Non-T"
seu$celltype2[seu$seurat_clusters %in% c(13)] <- "EEC"
seu$celltype2[seu$seurat_clusters %in% c(12)] <- "Pariet"
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "IM"
seu$celltype2 <- factor(seu$celltype2, levels = c("Non-T", "EEC", "Pariet", "IM", "Tumor"))
DimPlot(seu, group.by = "celltype2", cols = "alphabet", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes()
ggsave("celltype2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# save RDS ----
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
saveRDS(seu, file = "RDSfiles/seu_020_epi.RDS")

# UCell scoring ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")
seu <- DietSeurat(seu)
seu <- NormalizeData(seu)

H <- msigdbr(species    = "Homo sapiens", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list_human.RDS")
yap_list <- yap_list[4:5]

seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

## make ucell scores into z-scores
### change names for brevity
colnames(seu@meta.data)[colnames(seu@meta.data) == "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"] <- "EMT_UCell"
colnames(seu@meta.data)[colnames(seu@meta.data) == "INFLAMMATORY_RESPONSE_UCell"] <- "INFL_RESP_UCell"

ucell_cols <- grep("_UCell$", colnames(seu[[]]), value = TRUE)

for (cc in ucell_cols) {
  seu[[cc]] <- as.numeric(scale(seu@meta.data[[cc]]))
}

## create a composite score
seu[["TME_score"]] <- rowMeans(
  seu@meta.data[, c(
    "EMT_UCell",
    "HYPOXIA_UCell",
    "INFL_RESP_UCell"
  )],
  na.rm = TRUE
)



# FeatureScatter with YAP ----
## TME_score
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_TME.pdf", path = plot_path, width = 35, height = 30, units = "mm")

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_TME_no_p.pdf", path = plot_path, width = 35, height = 30, units = "mm", device = cairo_pdf)

## EMT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_EMT.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_EMT_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## HYPOXIA
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_HYPOXIA.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_HYPOXIA_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## INFL_RESP
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_INFL_RESP.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_INFL_RESP_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with G2M ----
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_G2M_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "dataset", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_G2M_TME_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# Deep vs Superficial ----
seu$layer <- factor(seu$layer, levels = c("S", "D"))
seu_cp <- seu
seu <- subset(seu_cp, subset = layer %in% c("S", "D"))
add_feat = "YAP_TARGETS_UCell"
VlnPlot(seu, group.by = "layer", features = add_feat, pt.size = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") &
  theme_panel() & NoLegend() &
  labs(title = "YAP_TARGETS", x = NULL, y = "z-score") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0(add_feat, "vln_", ".pdf"), path = plot_path, width = 25, height = 35, units = "mm")

for(pt in levels(seu$patient)){
  sub <- subset(seu, subset = seu$patient == pt)
  VlnPlot(sub, group.by = "layer", features = add_feat, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") &
    theme_panel() & NoLegend() &
    labs(title = "YAP_TARGETS", x = NULL, y = "z-score") &
    theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
  ggsave(paste0(add_feat, "vln_", pt, ".pdf"), path = plot_path, width = 25, height = 35, units = "mm")
}

add_feat = "TME_score"
VlnPlot(seu, group.by = "layer", features = add_feat, pt.size = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") &
  theme_panel() & NoLegend() &
  labs(title = "TME_score", x = NULL, y = "z-score") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0(add_feat, "vln_", ".pdf"), path = plot_path, width = 25, height = 35, units = "mm")

for(pt in levels(seu$patient)){
  sub <- subset(seu, subset = seu$patient == pt)
  VlnPlot(sub, group.by = "layer", features = add_feat, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") &
    theme_panel() & NoLegend() &
    labs(title = "TME_score", x = NULL, y = "z-score") &
    theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
  ggsave(paste0(add_feat, "vln_", pt, ".pdf"), path = plot_path, width = 25, height = 35, units = "mm")
}

