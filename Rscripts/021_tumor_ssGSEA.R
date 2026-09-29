# Continued from 020_
# 2026-03-20

# load packages ----
library(tidyverse)
library(openxlsx2)
library(Seurat)
library(msigdbr)
library(UCell)
library(ggpubr)
library(pals)
source("Rscripts/helpers.R")

# Make directories ----
analysis_step <- "021_tumor_ssGSEA"
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")
seu <- DietSeurat(seu)
seu <- subset(seu, subset = celltype2 %in% c("Tumor") & tissue_type != "NT")
seu <- NormalizeData(seu)


# UCell scoring ----
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
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_TME.pdf", path = plot_path, width = 35, height = 30, units = "mm")

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_TME_no_p.pdf", path = plot_path, width = 35, height = 30, units = "mm", device = cairo_pdf)

## EMT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_EMT.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_EMT_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## HYPOXIA
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_HYPOXIA.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_HYPOXIA_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## INFL_RESP
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_INFL_RESP.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_INFL_RESP_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with G2M ----
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_G2M_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_G2M_TME_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# Deep vs Superficial ----
seu$layer <- factor(seu$layer, levels = c("S", "D"))
levels(seu$layer) <- c("Superficial", "Deep")
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

