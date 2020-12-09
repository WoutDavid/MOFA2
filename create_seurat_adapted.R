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
counts <- Read10X_h5("/media/david/Puzzles/IBP/mouse/cellranger/mouse_filtered_feature_bc_matrix.h5")
#counts <- Read10X(paste0(basedir,"/e18_mouse_brain_fresh_5k_raw_feature_bc_matrix.h5"))

###################
## Create Seurat ##
###################
rm(seurat)
seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "scRNA+scATAC_mouse",
  min.cells = 1,
)

# Add ATAC modality, beware that this is still an AssayObject, not a chromatinAssay object
#this Assay Object is going to get transformed into two different chromatinAssay objects later on, in TrainingTheModel
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])

##################
## Add metadata ##
##################

#read in the metadata, replace the barcordes by barcodes without the -1 in there
#I decided not to replace the -1. It might be less pretty and interpretable, but it makes the merging
#of the metadata a lot easier
rm(metadata)
metadata <- fread("/media/david/Puzzles/IBP/mouse/cellranger/mouse_per_barcode_metrics.csv") %>%
 .[,barcode:=gsub("-1","-1",barcode)]

#filter those rows that do not have a 1 in their is_cell column
metadata <- metadata[metadata$is_cell==1]

dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  tibble::column_to_rownames("barcode")

##add metadata
seurat <- AddMetaData(seurat, dt)

###############################################################
## Quality control for mitochondrial genes and nfeatures_RNA ##
###############################################################
##watch out, the pattern for mice mitochondrial dna is lowercase mt, not uppercase, so this needs to be changed per dataset
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
VlnPlot(seurat, features = c("nFeature_RNA", "nFeature_ATAC", "percent.mt"), ncol = 3)
##for huma we choose: subset = nFeature_RNA < 10000 & percent.mt < 5 & nFeature_ATAC < 40000
#for mice we choos: subset(seurat, subset = nFeature_RNA < 7000 & percent.mt < 10 & nFeature_ATAC < 25000)
seurat <- subset(seurat, subset = nFeature_RNA < 7000 & percent.mt < 10 & nFeature_ATAC < 25000)


##########
## Save ##
##########
#optional ofcourse, takes some storage
saveRDS(seurat, "/media/david/Puzzles/IBP/mouse/third_model_mouse/third_create_seurat_mouse.RDS", compress = FALSE)
