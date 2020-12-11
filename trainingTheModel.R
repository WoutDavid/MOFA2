##################################################################################################################
## WATCH OUT: this script is reused for mouse/human seperately, so dont forget to replace all mentions of human ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################

#load and install all libraries:
library(data.table)
library(ggplot2)
library(Seurat)

#this makes it so that R is allowed to look on bioconducter to install packages
#it's allowing a lot more, but it almost makes sure your package gets installed
setRepositories(ind=1:2)
#didn't work, creator of Signac says you can also just install this manually:
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
#yup that worked

library(Signac)
library(msigdbr)
library(BiocManager)
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(TFBSTools)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

##choose the correct one for the correct organism:
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MOFA2)

#load Seurat object coming out of cell annotation
seurat_mofa_human <- readRDS("/media/david/Puzzles/IBP/final_model/Human_Annotated.Rds")
seurat_mofa_human

#downloading position specific weight matrix
#################################################################################################################################
##DONT FORGET TO CHANGE THIS TO THE CORRECT ORGANISM --> Mus musculus / Homo sapiens
##################################################################################################################################
pfm <- getMatrixSet(JASPAR2020,
                    opts = list(species = "Homo sapiens")
)

#import the feature metadata
feature_metadata <- fread("/media/david/Puzzles/IBP/human/cellranger/filtered_feature_bc_matrix/features.tsv") %>%
  setnames(c("ens_id","gene","view","chr","start","end"))

#extract rna metadatadata
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]

#extract atac metadatadata
#ens.null := null selects out that column, it just gets rid of it
feature_metadata.atac <- feature_metadata[view=="Peaks"] %>% 
  .[,ens_id:=NULL] %>% setnames("gene","peak")

#function to fix the naming problem in geneID
NameFix <- function(string){
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
##namefix is applied to the peak names of this document, because they don't match the original data, which makes merging impossible
foo$peak <- sapply(foo$peak, NameFix)

  
feature_metadata.atac <- feature_metadata.atac %>% 
  merge(foo,by="peak",all.x=TRUE)

#Split ATAC matrix depending on the peak type and create a ChromatinAssay for each modality using the Signac package.
#This object requires a GRanges object with the peak metadata. 
#BiocManager::install("GenomicRanges")
#BiocManager::install("motifmatchr")
library(motifmatchr)
library(GenomicRanges)
#first remove all peaks that are not on a chromosome
feature_metadata.atac <- feature_metadata.atac[startsWith(feature_metadata.atac$chr, "chr")==TRUE,]
#create ChromatinAssays
for (i in c("distal","promoter")) {
  # Create GRanges
  peaks.granges <- feature_metadata.atac %>%
    .[peak_type==i] %>%
    .[,c("chr","start","end","peak")] %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  # Scan motifs throughout the DNA sequence of each peak and create a binary matrix of motif-peak presence.
  ##DONT FORGET TO CHANGE THE GENOME OBJECT TO THE CORRECT ORGANISM
  motif.matrix <- CreateMotifMatrix(
    features = peaks.granges,
    pwm = pfm,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    use.counts = FALSE
  ) %>% as.matrix
  
  # AddChromatinAssay to the Seurat object
  seurat_mofa_human@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat_mofa_human@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    motifs = CreateMotifObject(motif.matrix, pfm)
  )
  
}


#normalizing and scaling
seurat_mofa_human <- NormalizeData(seurat_mofa_human, normalization.method = "LogNormalize", assay = "RNA")
seurat_mofa_human <- ScaleData(seurat_mofa_human, do.center = TRUE, do.scale = FALSE)
                                                                                                                                                                                                                                                                                                                                                        #normalizing the ATAC dataset
##we might want to try and regress out nFeature_ATAQ
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat_mofa_human <- RunTFIDF(seurat_mofa_human, assay = i)
}

#feature selection on the seurat objects to do preselection to make creating the mofa object easier
#the number of features we select here is definitely relevant, this might be worth looking into
seurat_mofa_human <- FindVariableFeatures(seurat_mofa_human, 
                               selection.method = "vst", 
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)

##FindtopFeatures with min.cutoff 1/5th of the cell population, it's the same ratio as in the vignette
##you ofcourse have to adapt this to the different datasets.
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat_mofa_human <- FindTopFeatures(seurat_mofa_human,min.cutoff = 1000, assay=i, verbose = FALSE)
  print(length(seurat_mofa_human[[i]]@var.features))
}

#save the seurat object before training the model, it can be an good source of extra information in the downstream analysis
saveRDS(seurat_mofa_human, "/media/david/Puzzles/IBP/final_model/human/seurat_human_final_model.RDS")

##Create mofa object
mofa <- create_mofa(seurat_mofa_human, assays = c("RNA","ATAC_distal","ATAC_promoter"))

##define the model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15

mofa <- prepare_mofa(mofa, model_options = model_opts)
##Careful: this command actually trains the model
mofa <- run_mofa(mofa, outfile = "/media/david/Puzzles/IBP/final_model/human/human_MOFA_final_model.hdf5", use_basilisk=TRUE)

##save mofa object for downstream analysis
saveRDS(mofa, "/media/david/Puzzles/IBP/final_model/human/mofa_object_human_final_model.RDS")
