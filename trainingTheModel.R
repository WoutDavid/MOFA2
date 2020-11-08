##################################################################################################################
## WATCH OUT: this script is reused for mouse/human seperately, so dont forget to replace all mentions of mouse ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################
#laod libraries:
library(data.table)
library(ggplot2)
library(Seurat)
#this makes it so that R is allowed to look on bioconducter to install packages
#it's allowing a lot more, but it almost makes sure your package gets installed
setRepositories(ind=1:2)
#didn't work, creator of Signac says you can also just install this:
# BiocManager::install(c(
#   'AnnotationFilter',
#   'BiocGenerics',
#   'GenomeInfoDb',
#   'GenomicFeatures',
#   'GenomicRanges',
#   'IRanges',
#   'Rsamtools',
#   'S4Vectors',
#   'TFBSTools',
#   'ggbio',
#   'motifmatchr',
#   'AnnotationDbi')
# )
#yup that worked, thank god

library(Signac)

library(msigdbr)
library(BiocManager)
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(TFBSTools)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Hsapiens.UCSC.hg38)

library(MOFA2)
seurat_mofa_mouse <- readRDS("/media/david/Puzzles/IBP/mouse/cellranger/seurat_mouse.RDS")
seurat_mofa_mouse
#doesnt work because metadata wasnt added
head(seurat_mofa_mouse@meta.data)

#Quality control
#todo: look at the suggestions send to us by the promotor

#downloading position specific weight matrix
pfm <- getMatrixSet(JASPAR2020,
                    opts = list(species = "Mus musculus")
)
pfm
#import the feature metadata
feature_metadata <- fread("/media/david/Puzzles/IBP/mouse/cellranger/filtered_feature_bc_matrix/features.tsv") %>%
  setnames(c("ens_id","gene","view","chr","start","end"))
feature_metadata
#extract rna data
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]
head(feature_metadata.rna)
tail(feature_metadata.rna)

#extract atac data
#ens.null := null selects out that column, it just gets rid of it
feature_metadata.atac <- feature_metadata[view=="Peaks"] %>% 
  .[,ens_id:=NULL] %>% setnames("gene","peak")
head(feature_metadata.atac)
tail(feature_metadata.atac)
#selfmade function to fix the naming problem in geneID
NameFix <- function(string){
  ##split the string by _
  words =strsplit(string, "_")[[1]]
  words <- append(words, ":", after = 1)
  words <- append(words, "-", after = 3)
  merged <- paste(words, collapse="")
  return(merged)
}

#import the atac annotations
foo <- fread("/media/david/Puzzles/IBP/mouse/cellranger/mouse_atac_peak_annotation.tsv") %>%
  .[,c("peak","peak_type")] %>%
  .[peak_type%in%c("distal", "promoter")]
head(foo)
tail(foo, 50)
foo$peak <- sapply(foo$peak, NameFix)
head(foo)
tail(foo)
##PROBLEM FOUND: column names dont match, some fumpbling with the names is in order
#should be fixed now, by using the function NameFix on the data coming from the atac_peak_annotation

  
feature_metadata.atac <- feature_metadata.atac %>% 
  merge(foo,by="peak",all.x=TRUE)
head(feature_metadata.atac)
tail(feature_metadata.atac)

#Split ATAC matrix depending on the peak type and create a ChromatinAssay for each modality using the Signac package.
#This object requires a GRanges object with the peak metadata. 
#this doesn't work yet, but it's supposed to divide the ATAC assay into a distal and a promoter atac assay
#problem is those last few x lines in the data that contain these sequence: "GL000194.1"
#We might wanna decide to just let them be
#BiocManager::install("GenomicRanges")
#BiocManager::install("motifmatchr")
library(motifmatchr)
library(GenomicRanges)

for (i in c("distal","promoter")) {
  # Create GRanges
  peaks.granges <- feature_metadata.atac %>%
    .[peak_type==i] %>%
    .[,c("chr","start","end","peak")] %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  # Scan motifs throughout the DNA sequence of each peak and create a binary matrix of motif-peak presence.
  # motif.matrix <- CreateMotifMatrix(
  #   features = peaks.granges,
  #   pwm = pfm,
  #   genome = BSgenome.Mmusculus.UCSC.mm10,
  #   use.counts = FALSE
  # ) %>% as.matrix
  
  # AddChromatinAssay to the Seurat object
  seurat_mofa_mouse@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat_mofa_mouse@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    #motifs = CreateMotifObject(motif.matrix, pfm)
  )
  
}
seurat_mofa_mouse

#alright so the 4 assays are there, but the motif matrix has not been added, so some downstream analysis will still not match
#let's prepare the data by normalizing and scaling
seurat_mofa_mouse <- NormalizeData(seurat_mofa_mouse, normalization.method = "LogNormalize", assay = "RNA")
seurat_mofa_mouse <- ScaleData(seurat_mofa_mouse, do.center = TRUE, do.scale = FALSE)

#normalizing the ATAC dataset
#some features contain 0 counts cause we didn't do any quality control yet
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat_mofa_mouse <- RunTFIDF(seurat_mofa_mouse, assay = i)
}

#feature selection on the seurat objects to do preselction to make creating the mofa thing easier
#the number of features we select here is definitly relevant, this might be worth looking into
seurat_mofa_mouse <- FindVariableFeatures(seurat_mofa_mouse, 
                               selection.method = "vst", 
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)
#this might be wrong, the vignette uses the commented out method FindTopFeatures, but that didn't select enough features out for the model to be trained
#seurat_mofa_mouse <- FindVariableFeatures(seurat_mofa_mouse, selection.method = "vst", nfeatures = 5000, assay="ATAC", verbose = FALSE)
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat_mofa_mouse <- FindVariableFeatures(seurat_mofa_mouse, assay=i, nfeatures = 5000, verbose = FALSE, selection.method = "vst")
  print(length(seurat_mofa_mouse[[i]]@var.features))
}
# for (i in c("ATAC_distal","ATAC_promoter")) {
#   seurat <- FindTopFeatures(seurat, assay=i, min.cutoff = 2000)
#   print(length(seurat[[i]]@var.features))
# }

#save the seurat object before training the model, it can be an good source of extra information in the downstream analysis
saveRDS(seurat_mofa_mouse, "seurat_mouse_before_training.RDS")
seurat_mofa_mouse
##Training the model woep woep!
mofa <- create_mofa(seurat_mofa_mouse, assays = c("RNA","ATAC_distal","ATAC_promoter"))
##again: notice here that there's only 5000 features in the ATAC seq datacause i used FindVariableFeatures, and not FindTopFeatures
mofa


model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
mofa
mofa <- prepare_mofa(mofa,
                     model_options = model_opts)
mofa <- run_mofa(mofa, outfile = "/media/david/Puzzles/IBP/mouse/mouse_MOFA_model.hdf5")
saveRDS(mofa, "/media/david/Puzzles/IBP/mouse/mofa_object_mouse.RDS")
