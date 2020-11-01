library(Seurat)
library(data.table)
library(purrr)
library(rhdf5)
basedir <- "/media/david/Puzzles/IBP/human/cellranger/filtered_feature_bc_matrix/"
outfile <- "/media/david/Puzzles/IBP/human/cellranger/seurat_human.RDS"


###############
## Load data ##
###############
counts <- Read10X_h5("/media/david/Puzzles/IBP/human/cellranger/human_filtered_feature_bc_matrix.h5")
#counts <- Read10X(paste0(basedir,"/e18_mouse_brain_fresh_5k_raw_feature_bc_matrix.h5"))

###################
## Create Seurat ##
###################

seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "scRNA+scATAC_human",
  min.cells = 1
)

# Add ATAC modality
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])

seurat

##################
## Add metadata ##
##################
##this isn't correct yet, since i don't have the correct metadata file
metadata <- fread("/media/david/Puzzles/IBP/human/cellranger/human_per_barcode_metrics.csv") %>%
  .[,barcode:=gsub("-1","",barcode)]

rename_celltypes <- c(
  "naive CD4 T cells" = "Lymphoid",
  "memory CD4 T cells" = "Lymphoid",
  "naive CD8 T cells" = "Lymphoid",
  "CD56 \\(bright\\) NK cells" = "Lymphoid",
  "CD56 \\(dim\\) NK cells" = "Lymphoid",
  "memory B cells" = "Lymphoid",
  "naive B cells" = "Lymphoid",
  "effector CD8 T cells" = "Lymphoid",
  "MAIT T cells" = "Lymphoid",
  "non-classical monocytes" = "Myeloid",
  "intermediate monocytes" = "Myeloid",
  "classical monocytes" = "Myeloid",
  "myeloid DC" = "Myeloid",
  "plasmacytoid DC" = "Lymphoid"
)
metadata$broad_celltype <- stringr::str_replace_all(metadata$celltype,rename_celltypes)

dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
  .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
  tibble::column_to_rownames("barcode")

seurat <- AddMetaData(seurat, metadata)

head(seurat@meta.data)

##########
## Save ##
##########
#note, you skipped a big part of the metadata cause the celltype thing didn't work 
saveRDS(seurat, outfile, compress = FALSE)
