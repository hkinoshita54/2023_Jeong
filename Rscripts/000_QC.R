# 2026-03-18

# Setting up ----
dataset <- "GSE167297"

## Make directories
analysis_step <- "000_QC"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(data.table)
library(Matrix)
library(Seurat)
source("Rscripts/helpers.R")

# load data ----
## SraRunTable downloaded from GEO's Sra run selector
sra_table <- read.delim(file = "data/SraRunTable.csv", sep = ",") 
meta <- sra_table %>% 
  select(gsm = Sample.Name, patient, tissue_type = tissue) %>% 
  distinct() %>% 
  arrange(gsm) %>% 
  mutate(
  layer = recode(tissue_type, "normal gastric tissue" = NA_character_, "superficial layer of gastric cancer" = "S", "deep layer of gastric cancer" = "D"),
  tissue_type = recode(tissue_type,"normal gastric tissue" = "NT","superficial layer of gastric cancer" = "T","deep layer of gastric cancer" = "T"),
  patient = sub("patient", "Jeong", patient)
)
nsample <- nrow(meta)

saveRDS(meta, "RDSfiles/metadata.RDS")

## Create a list of Seurat object
files_list <- list.files(path = "data/GSE167297_RAW", pattern = "\\.txt\\.gz$", full.names = T) %>% sort()

seu_list <- list()
for (f in files_list) {
  gsm <- sub("_.*$", "", basename(f))
  
  dt <- fread(f)
  genes <- make.unique(dt[[1]])
  mtx <- as.matrix(dt[, -1, with = FALSE])
  rownames(mtx) <- genes
  mtx <- Matrix(mtx, sparse = TRUE)
  
  sample_meta <- meta[meta$gsm == gsm, ]
  
  cell_meta <- data.frame(
    cell_ID = colnames(mtx),
    gsm = gsm
  ) %>%
    left_join(sample_meta, by = "gsm") %>%
    column_to_rownames("cell_ID")
  
  seu <- CreateSeuratObject(
    counts = mtx,
    project = gsm,
    min.cells = 3,
    min.features = 200,
    meta.data = cell_meta
  )
  
  seu_list[[gsm]] <- seu
}

## Sanity check
length(seu_list)
names(seu_list)

sapply(seu_list, ncol)
sapply(seu_list, function(x) unique(x$patient))
sapply(seu_list, function(x) unique(x$tissue_type))
sapply(seu_list, function(x) unique(x$layer))

## DoubletFinder
source("Rscripts/DoubletFinder.R")

## merge the list, clean metadata
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)], add.cell.ids = names(seu_list))
names(seu@meta.data)
seu@meta.data <- seu@meta.data[,c(1:7,9,11)]
str(seu@meta.data)

seu$dataset <- dataset
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu@meta.data)

# QC ----
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("QC_vln_unfil.pdf", path = plot_path, width = 75, height = 35, units = "mm")

seu <- subset(seu, subset = doublet_finder == "Singlet")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("QC_vln_filt.pdf", path = plot_path, width = 75, height = 35, units = "mm")

# Save filtered object
rna_layers <- Layers(seu[["RNA"]])
keep_layers <- rna_layers[grepl("^counts.", rna_layers)]
seu[["RNA"]] <- subset(seu[["RNA"]], layers = keep_layers)
seu <- JoinLayers(seu)
saveRDS(seu, file.path("RDSfiles", "seu_000_filt.RDS"))
