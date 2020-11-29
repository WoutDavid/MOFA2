#importing the seurat object and the MOFA model
library(Seurat)
library(MOFA2)
seurat <- readRDS("seurat_human.RDS")
mofa <- load_model("first_model.hdf5")

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
##Correlate factors & covariates
correlate_factors_with_covariates(mofa, 
                                  covariates = c("nFeature_RNA","nFeature_ATAC")
)

##Visualisation of factor values --> metadata needed
##plot_factor(mofa, factors=1, group_by = "celltype", color_by="broad_celltype") +
##  theme(
##    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
##  )

##Visualisation of feature weights
plot_weights(mofa, 
             view = "RNA", 
             factors = 1, 
             nfeatures = 20, 
             text_size = 4
)

##Visualisation covariation patterns in data (for factor 1) -->metadata
# plot_data_scatter(mofa, 
#                   view = "RNA", 
#                   factor = 1, 
#                   features = 6,
#                   color_by = "broad_celltype",
#                   add_lm = T,
#                   dot_size = 1
# )
# plot_top_weights(mofa, 
#                  view = "ATAC_distal", 
#                  factors = 1, 
#                  sign = "positive",
#                  nfeatures = 15,
# )
# plot_data_heatmap(mofa, 
#                   view = "ATAC_promoter", 
#                   factor = 1, 
#                   features = 50,
#                   show_rownames = F, show_colnames = F, 
#                   cluster_rows = T, cluster_cols = F,
#                   annotation_samples = "broad_celltype"
# )
# plot_data_heatmap(mofa, 
#                   view = "ATAC_promoter", 
#                   factor = 1, 
#                   features = 50,
#                   show_rownames = F, show_colnames = F, 
#                   cluster_rows = T, cluster_cols = F,
#                   annotation_samples = "broad_celltype",
#                   denoise = TRUE
# )

#Non-linear dimensionality reduction
##Using MOFA-factors
factors <- 1:get_dimensions(mofa)[["K"]]
factors <- factors[!factors%in%c(4,7)]

mofa <- run_umap(mofa, 
                 factors = factors, 
                 n_neighbors = 15,  
                 min_dist = 0.30
)

  #We need the celltype metadata to make it look nice
plot_dimred(mofa, 
            method = "UMAP", 
            color_by = "celltype", 
            label = TRUE, 
            stroke=0.05, 
            dot_size = 1, 
            legend = FALSE
) + scale_fill_manual(values=colors)

  #Visualising contribution of each factor for the different celltypes
for (i in paste0("Factor",1:3)) {
  p <- plot_dimred(mofa, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 1
  )
  print(p)
}


##Using only RNA-seq data
DefaultAssay(seurat) <- "RNA"

seurat_mofa_human <- FindVariableFeatures(seurat,
                                          selection.method = "vst",
                                          nfeatures = 5000,
                                          assay = "RNA",
                                          verbose = FALSE
)
seurat_mofa_human <- NormalizeData(seurat_mofa_human, normalization.method = "LogNormalize", assay = "RNA")
seurat_mofa_human <- ScaleData(seurat_mofa_human, do.center = TRUE, do.scale = FALSE)


seurat <- RunPCA(seurat_mofa_human, npcs = 50, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'pca', dims = 1:50, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes()


##Using only ATAC
DefaultAssay(seurat) <- "ATAC_distal"
seurat <- RunSVD(seurat, n = K, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:K, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
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








