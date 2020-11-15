##################################################################################################################
## WATCH OUT: this script is reused for mouse/human seperately, so dont forget to replace all mentions of mouse ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################

library(Seurat)
library(data.table)
library(purrr)
library(rhdf5)


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
  min.cells = 1,
)
seurat
dim(seurat@assays$RNA@counts)

# Add ATAC modality, beware that this is still an AssayObject, not a chromatinAssay object
#this Assay Object is going to get transformed into two different chromatinAssay objects later on, in TrainingTheModel
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])
dim(seurat@assays$ATAC@counts)
##86394 features, 3332 cells
seurat



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

dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  tibble::column_to_rownames("barcode")
dt
seurat <- AddMetaData(seurat, dt)
head(seurat@meta.data)
dim(seurat@meta.data)


###############################################################
## Quality control for mitochondrial genes and nfeatures_RNA ##
###############################################################
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##10 is not yet a very strict threshold for mitochondrial percentage, neither is 13000 feature_rna
seurat <- subset(seurat, subset = nFeature_RNA < 13000 & percent.mt < 10)

##########
## Save ##
##########
saveRDS(seurat, "second_create_seurat.RDS", compress = FALSE)
