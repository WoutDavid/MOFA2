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
rm(seurat)
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

rm(metadata)
test <-fread("/media/david/Puzzles/IBP/human/cellranger/human_per_barcode_metrics.csv")
dim(test)
#read in the metadata, replace the barcordes by barcodes without the -1 in there
metadata <- fread("/media/david/Puzzles/IBP/human/cellranger/human_per_barcode_metrics.csv") %>%
 .[,barcode:=gsub("-1","",barcode)]
dim(metadata)
#filter those rows that do not have a 1 in their is_cell column
metadata <- metadata[metadata$is_cell==1]
dim(metadata)
#this is code that only works if you have this specific metadata, it's not necessary.
# dt <- data.table(barcode=colnames(seurat)) %>%
#   merge(metadata,by="barcode", all.x=TRUE) %>%
#   .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
#   .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
#   tibble::column_to_rownames("barcode")


seurat <- AddMetaData(seurat, metadata)
?AddMetaData
head(seurat@meta.data)
dim(seurat@meta.data)
sum(is.na(seurat@meta.data$is_cell))
##########
## Save ##
##########
#note, you skipped a big part of the metadata cause the celltype thing didn't work 
saveRDS(seurat, outfile, compress = FALSE)
