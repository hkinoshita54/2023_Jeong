library(tidyverse)
library(Seurat)

harmonize_symbols <- function(seu, assay = "RNA", species = "human"){
  counts <- LayerData(seu, assay = assay, layer = "counts")
  chk <- HGNChelper::checkGeneSymbols(rownames(counts), species = species, unmapped.as.na = FALSE)
  new_symbols <- chk$Suggested.Symbol
  
  f <- factor(new_symbols, levels = unique(new_symbols))
  G <- Matrix::sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(rownames(counts), levels(f))
  )
  collapsed <- Matrix::t(G) %*% counts
  
  seu <- CreateSeuratObject(counts = collapsed, meta.data = seu@meta.data)
  
  return(seu)
}
