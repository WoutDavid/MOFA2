library(Seurat)
library(data.table)
library(purrr)
library(rhdf5)
basedir <- "/media/david/Puzzles/IBP/mouse/cellranger/"
outfile <- "/media/david/Puzzles/IBP/mouse/cellranger/seurat_mouse.RDS"


###############
## Load data ##
###############
rm(counts)
counts <- Read10X_h5("/media/david/Puzzles/IBP/mouse/cellranger/mouse_filtered_feature_bc_matrix.h5")
#counts <- Read10X(paste0(basedir,"/e18_mouse_brain_fresh_5k_raw_feature_bc_matrix.h5"))

###################
## Create Seurat ##
###################
rm(seurat)
seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "scRNA+scATAC_mouse",
  min.cells = 1
)
seurat
seurat@assays$RNA@counts
rownames(seurat@assays$RNA@counts)
colnames(seurat@assays$RNA@counts)
dim(seurat@assays$RNA@counts)

# Add ATAC modality
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])
dim(seurat@assays$ATAC@counts)
rownames(seurat@assays$ATAC@counts)
colnames(seurat@assays$ATAC@counts)

##all gene names are in caps
##all cells are in caps with regex: CAAGACAAGGACCTTG-1
seurat

##################
## Add metadata ##
##################

rm(metadata)
test <-fread("/media/david/Puzzles/IBP/human/cellranger/human_per_barcode_metrics.csv")
dim(test)
colnames(test)
rownames(test)
#read in the metadata, replace the barcordes by barcodes without the -1 in there
metadata <- fread("/media/david/Puzzles/IBP/mouse/cellranger/mouse_per_barcode_metrics.csv") %>%
  .[,barcode:=gsub("-1","-1",barcode)]
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
dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  tibble::column_to_rownames("barcode")

seurat <- AddMetaData(seurat, dt)
head(seurat@meta.data)
dim(seurat@meta.data)
sum(is.na(seurat@meta.data$is_cell))
##########
## Save ##
##########
#note, you skipped a big part of the metadata cause the celltype thing didn't work 
saveRDS(seurat, outfile, compress = FALSE)
