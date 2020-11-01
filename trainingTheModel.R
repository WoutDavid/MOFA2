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
head(seurat_mofa_human@meta.data[,c("celltype","broad_celltype","pass_rnaQC","pass_accQC")])
#same
table(seurat_mofa_human@meta.data$celltype)
#same
table(seurat_mofa_human@meta.data$broad_celltype)
#Quality control
seurat_PBMC <- seurat_PBMC %>%
  .[,seurat_PBMC@meta.data$pass_accQC==TRUE & seurat_PBMC@meta.data$pass_rnaQC==TRUE]

#downloading position specific weight matrix
pfm <- getMatrixSet(JASPAR2020,
                    opts = list(species = "Homo sapiens")
)

#import the metadata
feature_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/filtered_feature_bc_matrix/features.tsv.gz") %>%
  setnames(c("ens_id","gene","view","chr","start","end"))

#extract rna data
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]
#extract atac data
feature_metadata.atac <- feature_metadata[view=="Peaks"] %>% 
  .[,ens_id:=NULL] %>% setnames("gene","peak")
#import the atac annotations
foo <- fread("/media/david/Puzzles/IBP/PBMC/atac_peak_annotation.tsv") %>%
  .[,c("peak","peak_type")] %>%
  .[peak_type%in%c("distal", "promoter")]

feature_metadata.atac <- feature_metadata.atac %>% 
  merge(foo,by="peak",all.x=TRUE)

#Split ATAC matrix depending on the peak type and create a ChromatinAssay for each modality using the Signac package.
#This object requires a GRanges object with the peak metadata. 
#this doesnt work yet, but it's supposed to divide the ATAC assay into a distal and a promotor atac assay
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
  seurat_PBMC@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat_PBMC@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    motifs = CreateMotifObject(motif.matrix, pfm)
  )
  
}
seurat_PBMC

#alright so normally we should have 4 assays by now, the split didn't work but i'm gonna try it with two
#let's prepare the data by normalizing and scaling
seurat_PBMC <- NormalizeData(seurat_PBMC, normalization.method = "LogNormalize", assay = "RNA")
seurat_PBMC <- ScaleData(seurat_PBMC, do.center = TRUE, do.scale = FALSE)

seurat_PBMC <- RunTFIDF(seurat_PBMC, assay = "ATAC")

#feature selection on the seurat objects to do preselction to make creating the mofa thing easier?
seurat_PBMC <- FindVariableFeatures(seurat_PBMC, 
                               selection.method = "vst", 
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)

seurat_PBMC <- FindTopFeatures(seurat_PBMC, assay="ATAC", min.cutoff = 2000)

##Training the model woep woep!
mofa <- create_mofa(seurat_PBMC, assays = c("RNA","ATAC"))
mofa

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15

mofa <- prepare_mofa(mofa,
                     model_options = model_opts)
