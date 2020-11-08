##################################################################################################################
## WATCH OUT: this script is reused for mouse/human seperately, so dont forget to replace all mentions of mouse ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################

library(Seurat)
library(data.table)
library(purrr)
library(rhdf5)
basedir <- "/media/david/Puzzles/IBP/human/cellranger/filtered_feature_bc_matrix/"
outfile <- "/media/david/Puzzles/IBP/human/cellranger/seurat_human.RDS"


###############
## Load data ##
###############
rm(counts)
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
seurat
seurat@assays$RNA@counts
dim(seurat@assays$RNA@counts)

# Add ATAC modality
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])
dim(seurat@assays$ATAC@counts)
##86394 features, 3332 cells

##all gene names are in caps
##all cells are in caps with regex: CAAGACAAGGACCTTG-1
seurat
#rows are cells
##################
## Add metadata ##
##################



#read in the metadata, replace the barcordes by barcodes without the -1 in there
#I decided not to replace the -1. It might be less pretty and interpretable, but it makes the merging
#of the metadata a lot easier
rm(metadata)
metadata <- fread("/media/david/Puzzles/IBP/human/cellranger/human_per_barcode_metrics.csv") %>%
 .[,barcode:=gsub("-1","-1",barcode)]
dim(metadata)
head(metadata)
#filter those rows that do not have a 1 in their is_cell column
metadata <- metadata[metadata$is_cell==1]
dim(metadata)

#this is code that only works if you have this specific metadata, it's not necessary.
# dt <- data.table(barcode=colnames(seurat)) %>%
#   merge(metadata,by="barcode", all.x=TRUE) %>%
#   .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
#   .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
#   tibble::column_to_rownames("barcode")

head(metadata[,1])
head(seurat@meta.data)
dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  tibble::column_to_rownames("barcode")
dt
seurat <- AddMetaData(seurat, dt)
head(seurat@meta.data)
dim(seurat@meta.data)

##########
## Save ##
##########
saveRDS(seurat, outfile, compress = FALSE)
