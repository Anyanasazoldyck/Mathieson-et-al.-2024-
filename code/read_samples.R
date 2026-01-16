#upload data as a Seurat object ####




library(Seurat)
library(stringr)
library(optparse)


## ---- Parse arguments ----
option_list <- list(
  make_option("--data_dir", type = "character"),
  make_option("--out_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

files <- list.files(
  opt$data_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

## ---- Initialize list ----
sc_list <- vector("list", length(files))
names(sc_list) <- basename(files)

## ---- Loop ----
for (i in seq_along(files)) {
  
  sample <- files[i]
  sample_id <- str_extract(basename(sample), "P[0-9]+")
  
  expr <- read.table(
    sample,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(expr),
    assay = "RNA",
    project = sample_id
  )
  
  seurat_obj$patient_id <- sample_id
  sc_list[[i]] <- seurat_obj
}

## ---- Save ----
saveRDS(
  sc_list,
  file = file.path(opt$out_dir, "late_nsclc_raw_list.rds")
)
