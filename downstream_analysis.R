#importing the seurat object and the MOFA model
library(Seurat)
library(MOFA2)
seurat <- readRDS("seurat_human_final_model.RDS")
mofa <- load_model("human_MOFA_final_model.hdf5")

#Add metadata
library(data.table)
samples_metadata(mofa) <- seurat@meta.data %>%
  tibble::rownames_to_column("sample") %>%
  as.data.table

#Assess correlation
plot_factor_cor(mofa)

#Variance decomposition
plot_variance_explained(mofa, max_r2 = 4)
plot_variance_explained(mofa, plot_total = TRUE)[[2]]

#Characterization of factors
install.packages("psych")
library(psych)
library(ggplot2)
##Correlate factors & covariates
correlate_factors_with_covariates(mofa, 
                                  covariates = c("nFeature_RNA","nFeature_ATAC")
)

#Visualisation of factor values --> metadata needed
plot_factor(mofa, factors=1, group_by = "cell_type", color_by="cell_type") +
 theme(
   axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
 )


##Visualisation of feature weights
plot_weights(mofa, 
             view = "RNA", 
             factors = 1, 
             nfeatures = 20, 
             text_size = 4
)

##Visualisation covariation patterns in data (for factor 1)
install.packages("ggpubr")
library(ggpubr)
plot_data_scatter(mofa,
                  view = "RNA",
                  factor = 1,
                  features = 6,
                  color_by = "cell_type",
                  add_lm = T,
                  dot_size = 2
)
plot_top_weights(mofa,
                 view = "ATAC_distal",
                 factors = 1,
                 sign = "positive",
                 nfeatures = 15,
)
plot_data_heatmap(mofa,
                  view = "ATAC_promoter",
                  factor = 1,
                  features = 50,
                  show_rownames = F, show_colnames = F,
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "cell_type"
)

plot_data_heatmap(mofa,
                  view = "ATAC_promoter",
                  factor = 1,
                  features = 50,
                  show_rownames = F, show_colnames = F,
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "cell_type",
                  denoise = TRUE
)

#Non-linear dimensionality reduction
##Using MOFA-factors
factors <- 1:get_dimensions(mofa)[["K"]]
factors <- factors[!factors%in%c(4,7)]

mofa <- run_umap(mofa, 
                 factors = factors, 
                 n_neighbors = 15,  
                 min_dist = 0.30
)

  
colors<- c("red", "orange", "yellow", "green", "lightblue", "darkblue", "purple")
plot_dimred(mofa, 
            method = "UMAP", 
            color_by = "cell_type", 
            label = TRUE, 
            stroke=0, 
            dot_size = 0.5, 
            legend = FALSE
) + scale_fill_manual(values=colors)

  #Visualising contribution of each factor for the different celltypes
for (i in paste0("Factor",1:3)) {
  p <- plot_dimred(mofa, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 3
  )
  print(p)
}


##Using only RNA-seq data
DefaultAssay(seurat) <- "RNA"

#seurat_mofa_human <- FindVariableFeatures(seurat,
#                                          selection.method = "vst",
#                                          nfeatures = 5000,
#                                          assay = "RNA",
#                                          verbose = FALSE
#)
#seurat_mofa_human <- NormalizeData(seurat_mofa_human, normalization.method = "LogNormalize", assay = "RNA")
#seurat_mofa_human <- ScaleData(seurat_mofa_human, do.center = TRUE, do.scale = FALSE)

ElbowPlot(seurat, reduction='pca')
##Run with 6 PCs
seurat <- RunPCA(seurat_mofa_human, npcs = 6, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'pca', dims = 1:6, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes()+ scale_fill_manual(values=colors)

##Using only ATAC
DefaultAssay(seurat) <- "ATAC_distal"
seurat <- RunSVD(seurat, n = 6, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:6, verbose = FALSE)
DimPlot(seurat, label = T, reduction="umap") + 
  NoLegend() + NoAxes() + scale_fill_manual(values=colors)



##############
#GSEA mouse ##
##############

##install packages for the Gene sets
library(msigdbr)
library(cowplot)
#BiocManager::install("MOFAdata")
library(MOFAdata)

##load data
data("MSigDB_v6.0_C5_mouse")
head(colnames(MSigDB_v6.0_C5_mouse))


############## dont forget to change to the correct feature set

##To match the gene names,wwe have to capitalize the features
features_names(mofa_mouse)[["RNA"]] <- toupper(features_names(mofa_mouse)[["RNA"]])
features_names(mofa_mouse)[["RNA"]]
gsea.positive <- run_enrichment(mofa_mouse, 
                                feature.sets = MSigDB_v6.0_C5_mouse, 
                                view = "RNA",
                                sign = "positive"
)
# GSEA on negative weights
gsea.negative <- run_enrichment(mofa_mouse, 
                                feature.sets = MSigDB_v6.0_C5_mouse, 
                                view = "RNA",
                                sign = "negative"
)
names(gsea.positive)

##visualisation
##POSITIVE
#check if all factors have some enriched pathways:
plot_enrichment_heatmap(gsea.positive)

plot_enrichment(gsea.positive, factor = 1, max.pathways = 15)

plot_enrichment_detailed(gsea.positive,
                         factor = 1,
                         max.genes = 10,
                         max.pathways = 5
)

##we see some interesting genes, that you can fill in this list and then visualise their GSEA specifically:
genes <- list("x","x")

genes %>% map(~ plot_factors(mofa_mouse, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

##NEGATIVE
plot_enrichment_heatmap(gsea.negative)

plot_enrichment(gsea.negative, factor = 1, max.pathways = 15)

plot_enrichment_detailed(gsea.negative,
                         factor = 1,
                         max.genes = 10,
                         max.pathways = 5
)
##we see some interesting genes again in the negtive
genes <- list("x","x")

genes %>% map(~ plot_factors(mofa_mouse, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)



##############
#GSEA human ##
##############

data("MSigDB_v6.0_C5_human")

#following is commented out because it requires manual steps that I already performed, you can continue to the part where
##you read in the gene_names
# identifiers <- colnames(MSigDB_v6.0_C5_human)
# lapply(identifiers, write, "test.txt", append=TRUE)

## the result of the python code is Gene_names.txt, which is too large to put in the github, so you'll have to download
#it from google drive. I'll put it in the second_model directory.

#We have to integrate that back into the matrix:
gene_names <-  read.table("Gene_names.txt", sep="\t")
##note that there are empty entries in there, not all features are recognised but that shouldn't be a problem
vector <- gene_names$V2
colnames(MSigDB_v6.0_C5_human) <- vector

##To match the gene names,wwe have to capitalize the features
features_names(mofa)[["RNA"]] <- toupper(features_names(mofa)[["RNA"]])
features_names(mofa)[["RNA"]]
gsea.positive <- run_enrichment(mofa, 
                                feature.sets = MSigDB_v6.0_C5_human, 
                                view = "RNA",
                                sign = "positive"
)
# GSEA on negative weights
gsea.negative <- run_enrichment(mofa, 
                                feature.sets = MSigDB_v6.0_C5_human, 
                                view = "RNA",
                                sign = "negative"
)
names(gsea.positive)

##visualisation
##POSITIVE
#check if all factors have some enriched pathways:
plot_enrichment_heatmap(gsea.positive)

plot_enrichment(gsea.positive, factor = 1, max.pathways = 15)

plot_enrichment_detailed(gsea.positive,
                         factor = 1,
                         max.genes = 10,
                         max.pathways = 5
)

##we see some interesting genes, that you can fill in this list and then visualise their GSEA specifically:
genes <- list("x","x")

genes %>% map(~ plot_factors(mofa_human, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

##NEGATIVE
plot_enrichment_heatmap(gsea.negative)

plot_enrichment(gsea.negative, factor = 1, max.pathways = 15)

plot_enrichment_detailed(gsea.negative,
                         factor = 1,
                         max.genes = 10,
                         max.pathways = 5
)
##we see some interesting genes again in the negtive
genes <- list("x","x")

genes %>% map(~ plot_factors(mofa_human, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

##########################
#Motif enrichment human ##
##########################
##we perform motif enrichment on distal ATAC, because it explains the most variance
library(Seurat)
library(MOFA2)
library(Signac)

motif.matrix <- t(as.matrix(seurat[["ATAC_distal"]]@motifs@data))
rownames(motif.matrix)
motif.enrichment.positive <- run_enrichment(mofa,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "positive"
)

motif.enrichment.negative <- run_enrichment(mofa_human,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "negative"
)

##Visualising
###this one doesn't say much cause the motifs have bad interpretable names
plot_enrichment(motif.enrichment.positive, factor = 1, max.pathways = 15)

plot_enrichment(motif.enrichment.negative, factor = 1, max.pathways = 15)

##this is pretty cool, it uses the Signac package
sig.motifs.positive <- motif.enrichment.positive$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat[["ATAC_distal"]], motifs = sig.motifs.positive)

sig.motifs.negative <- motif.enrichment.negative$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat_human[["ATAC_distal"]], motifs = sig.motifs.negative)

############################
## Motif enrichment mouse ##
############################

motif.matrix <- t(as.matrix(seurat_mouse[["ATAC_distal"]]@motifs@data))
#####debug
motif.enrichment.positive <- run_enrichment(mofa_mouse,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "positive"
)

motif.enrichment.negative <- run_enrichment(mofa_mouse,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "negative"
)

##Visualising
###this one doesn't say much cause the motifs have bad interpretable names
plot_enrichment(motif.enrichment.positive, factor = 1, max.pathways = 15)

plot_enrichment(motif.enrichment.negative, factor = 1, max.pathways = 15)

##this is pretty cool, it uses the Signac package
sig.motifs.positive <- motif.enrichment.positive$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat_mouse[["ATAC_distal"]], motifs = sig.motifs.positive)

sig.motifs.negative <- motif.enrichment.negative$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat_mouse[["ATAC_distal"]], motifs = sig.motifs.negative)

