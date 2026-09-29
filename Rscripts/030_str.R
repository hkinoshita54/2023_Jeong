# Continued from 010_
# 2026-03-20
# 2026-05-26

## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "030_str"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))
library(openxlsx2)
library(msigdbr)
library(UCell)
library(ggpubr)
library(pals)

# Helper functions ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 6) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering ----
# Clustering w/o integration ----
seu <- readRDS("RDSfiles/seu_010_all.RDS")
seu <- subset(seu, subset = celltype1 %in% c("Fibro", "EC"))

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
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.pdf", path = plot_path, width = 35, height = 30, units = "mm")

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

add_feat <- c("COL1A1", "DCN", "ACTA2", "MYH11", "RGS5", "VWF", "LYVE1", "PTPRC", "EPCAM")
sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

## removing immune contamination, then cluster again
seu <- subset(seu, subset = seurat_clusters %in% c(8,10,13,16,18:22), invert = TRUE)

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient_ID)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Cluster") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "patient_ID", raster = TRUE, pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("pt2.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "doublet_finder", raster = TRUE, pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "doublet")
ggsave("pt2.pdf", path = plot_path, width = 25, height = 30, units = "mm")


seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers2.csv"))

## annotation
seu <- subset(seu, subset = seurat_clusters %in% c(12), invert = TRUE) # likely contamination
seu$celltype2 <- "Fibro"
seu$celltype2[seu$seurat_clusters %in% c(0,1,2,6,7,9,13)] <- "EC"
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "Myo"
seu$celltype2 <- factor(seu$celltype2, levels = c("Fibro", "EC", "Myo"))
DimPlot(seu, group.by = "celltype2", cols = "alphabet", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes()
ggsave("celltype2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue3.pdf", path = plot_path, width = 35, height = 30, units = "mm")

## dotplot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))
features <- c(
  "PDGFRA", 
  "LUM", 
  "PLVAP", 
  "VWF",
  "RGS5", 
  "ACTA2"
)
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  )+
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 35, height = 35, units = "mm")

# save RDS ----
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu$RNA_snn_res.1 <- NULL
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
  celltype2,
  seurat_clusters
)
saveRDS(seu, file = "RDSfiles/seu_030_str.RDS")
