# Continued from 010_
# 2026-06-02

## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "060_b_plasma"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))
library(openxlsx2)

# Helper functions ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 6) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering ----
seu <- readRDS("RDSfiles/seu_010_all.RDS")
seu <- subset(seu, subset = celltype1 %in% c("Plasma", "Bcell"))

## w/ harmony
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000)
hvg <- VariableFeatures(seu)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 30
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster_res0.5.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# Check markers ----
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.5.csv"))

## removing contamination
seu <- subset(seu, subset = seurat_clusters %in% c(7), invert = TRUE) # obvious contamination

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster2.pdf", path = plot_path, width = 35, height = 30, units = "mm")


# save RDS ----
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu$sample <- seu$orig.ident
seu@meta.data <- seu@meta.data %>% select(
  orig.ident,
  dataset,
  nCount_RNA,
  nFeature_RNA,
  sample,
  patient,
  tissue_type,
  layer,
  percent.mt,
  S.Score,
  G2M.Score,
  Phase,
  celltype1,
  seurat_clusters
)
saveRDS(seu, file = "RDSfiles/seu_060_b_plasma.RDS")
