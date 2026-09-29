# copykat
# 2026-03-27

# Setting up ----

## Make directories
analysis_step <- "100_copykat"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
library(copykat)
source("Rscripts/helpers.R")
source("Rscripts/harmonize_symbols.R")

# load data ----
epi <- readRDS("RDSfiles/seu_020_epi.RDS")
str <- readRDS("RDSfiles/seu_030_str.RDS")
# mye <- readRDS("RDSfiles/seu_040_mye.RDS")
seu <- merge(x = epi, y = str)
# seu <- merge(x = epi, y = list(str, mye))
seu <- JoinLayers(seu)
seu <- harmonize_symbols(seu)
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
seu_all <- seu
seu_tum <- subset(seu_all, subset = tissue_type == "T")
table(seu_tum$patient, seu_tum$celltype1)

# copykat per patient ----
seu_list <- SplitObject(seu_tum, split.by = "patient")
for(pt in names(seu_list)){
  message("Processing: ", pt)
  
  tryCatch({
    seu <- seu_list[[pt]]
    mtx <- LayerData(seu, assay = "RNA", layer = "counts")
    dip_cells <- colnames(seu)[seu$celltype1 != "Epi"]
    
    set.seed(1)
    ck <- withr::with_dir(
      res_path,
      copykat(
        rawmat = mtx,
        id.type = "S",
        ngene.chr = 5,
        win.size = 25,
        KS.cut = 0.10,
        sam.name = pt,
        distance = "euclidean",
        norm.cell.names = dip_cells,
        output.seg = FALSE,
        plot.genes = FALSE,
        genome = "hg20",
        n.cores = 8
      )
    )
    
    pred <- as.data.frame(ck$prediction)
    seu$copykat <- ""
    seu$copykat[match(pred$cell.names, colnames(seu))] <- pred$copykat.pred
    seu_list[[pt]] <- seu
    
    seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
    npcs <- 20
    seu <- RunPCA(seu, npcs = npcs)
    seu <- FindNeighbors(seu, dims = 1:npcs)
    seu <- FindClusters(seu, resolution = 0.5)
    n_neighbors <- min(30, ncol(seu) - 1)
    seu <- RunUMAP(seu, dims = 1:npcs, n.neighbors = n_neighbors)
    
    DimPlot(seu, group.by = "celltype1", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
      theme_panel() & NoAxes() & labs(title = "celltype1") 
    ggsave(paste0(pt, "_celltype1.pdf"), path = plot_path, width = 30, height = 30, units = "mm")
    
    DimPlot(seu, group.by = "copykat", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
      theme_panel() & NoAxes() & labs(title = "copykat")
    ggsave(paste0(pt, "_copykat.pdf"), path = plot_path, width = 35, height = 30, units = "mm")
  }, error = function(e) {
    message("ERROR in ", pt, ": ", e$message)
  })
}

seu_tum <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
seu_tum <- JoinLayers(seu_tum)
table(seu_tum$celltype1, seu_tum$copykat)
copykat_res <- seu_tum[[]]["copykat"]
saveRDS(copykat_res, "RDSfiles/copykat_res_ref_str.RDS")

tum <- subset(seu_tum, subset = celltype1 == "Epi")
saveRDS(tum, "RDSfiles/seu_100_tum.RDS")
