#laod libraries:
library(data.table)
library(ggplot2)
library(Seurat)
#this makes it so that R is allowed to look on bioconducter to install packages
#it's allowing a lot more, but it almost makes sure your package gets installed
setRepositories(ind=1:2)
#didn't work, creator of Signac says you can also just install this:
BiocManager::install(c(
  'AnnotationFilter',
  'BiocGenerics',
  'GenomeInfoDb',
  'GenomicFeatures',
  'GenomicRanges',
  'IRanges',
  'Rsamtools',
  'S4Vectors',
  'TFBSTools',
  'ggbio',
  'motifmatchr',
  'AnnotationDbi')
)
#yup that worked, thank god
library(Signac)

library(msigdbr)
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(TFBSTools)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

library(MOFA2)
seurat_mofa_human <- readRDS("/media/david/Puzzles/IBP/human/cellranger/seurat_human.RDS")
seurat_mofa_human
#doesnt work because metadata wasnt added
head(seurat_mofa_human@meta.data)

#Quality control
#todo: look at the suggestions send to us by the promotor

#downloading position specific weight matrix
pfm <- getMatrixSet(JASPAR2020,
                    opts = list(species = "Homo sapiens")
)
pfm
#import the feature metadata
feature_metadata <- fread("/media/david/Puzzles/IBP/human/cellranger/filtered_feature_bc_matrix/features.tsv") %>%
  setnames(c("ens_id","gene","view","chr","start","end"))

#extract rna data
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]
#extract atac data
#ens.null := null selects out that column, it just gets rid of it
feature_metadata.atac <- feature_metadata[view=="Peaks"] %>% 
  .[,ens_id:=NULL] %>% setnames("gene","peak")

#selfmade function to fix the naming problem in geneID
nameFix <- function(string){
  ##split the string by _
  words =strsplit(string, "_")[[1]]
  words <- append(words, ":", after = 1)
  words <- append(words, "-", after = 3)
  merged <- paste(words, collapse="")
  return(merged)
}

#import the atac annotations
foo <- fread("/media/david/Puzzles/IBP/human/cellranger/human_atac_peak_annotation.tsv") %>%
  .[,c("peak","peak_type")] %>%
  .[peak_type%in%c("distal", "promoter")]
foo$peak
foo$peak <- sapply(foo$peak, nameFix)
foo$peak
##PROBLEM FOUND: column names dont match, some fumpbling with the names is in order
#should be fixed now, except for the last 10 or so entries that don't follow the same regime, i'm ignoring those for now

  
feature_metadata.atac <- feature_metadata.atac %>% 
  merge(foo,by="peak",all.x=TRUE)
feaute_metadata.atac
#Split ATAC matrix depending on the peak type and create a ChromatinAssay for each modality using the Signac package.
#This object requires a GRanges object with the peak metadata. 
#this doesnt work yet, but it's supposed to divide the ATAC assay into a distal and a promotor atac assay
#problem is those last few x lines in the data that contain these sequence: "GL000194.1"
#We might wanna decide to just let them be
BiocManager::install("GenomicRanges")
BiocManager::install("motifmatchr")
library(motifmatchr)
library(GenomicRanges)
for (i in c("distal","promoter")) {
  # Create GRanges
  peaks.granges <- feature_metadata.atac %>%
    .[peak_type==i] %>%
    .[,c("chr","start","end","peak")] %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  # Scan motifs throughout the DNA sequence of each peak and create a binary matrix of motif-peak presence.
  motif.matrix <- CreateMotifMatrix(
    features = peaks.granges,
    pwm = pfm,
    genome = 'hg38',
    use.counts = FALSE
  ) %>% as.matrix
  
  # AddChromatinAssay to the Seurat object
  seurat_mofa_human@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat_mofa_human@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    motifs = CreateMotifObject(motif.matrix, pfm)
  )
  
}
seurat_mofa_human

#alright so normally we should have 4 assays by now, the split didn't work but i'm gonna try it with two
#let's prepare the data by normalizing and scaling
seurat_mofa_human <- NormalizeData(seurat_mofa_human, normalization.method = "LogNormalize", assay = "RNA")
seurat_mofa_human <- ScaleData(seurat_mofa_human, do.center = TRUE, do.scale = FALSE)

#normalizing the ATAC dataset
#soem features contain 0 counts cause we didn't do any quality control yet
seurat_mofa_human <- RunTFIDF(seurat_mofa_human, assay = "ATAC")

#feature selection on the seurat objects to do preselction to make creating the mofa thing easier
#the number of features we select here is definitly relevant, this might be worth looking into
seurat_mofa_human <- FindVariableFeatures(seurat_mofa_human, 
                               selection.method = "vst", 
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)
#this might be wrong, the vignette uses the commented out method FindTopFeatures, but that didn't select enough features out for the model to be trained
seurat_mofa_human <- FindVariableFeatures(seurat_mofa_human, selection.method = "vst", nfeatures = 5000, assay="ATAC", verbose = FALSE)
#seurat_mofa_human <- FindTopFeatures(seurat_mofa_human, assay="ATAC", min.cutoff = 2000)

#save the seurat object before training the model, it can be an good source of extra information in the downstream analysis
saveRDS(seurat_mofa_human, "seurat_human_before_training.RDS")
##Training the model woep woep!
library(MOFA2)
mofa <- create_mofa(seurat_mofa_human, assays = c("RNA","ATAC"))


model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
mofa
mofa <- prepare_mofa(mofa,
                     model_options = model_opts)
mofa <- run_mofa(mofa, outfile = "MOFA_model.hdf5")
