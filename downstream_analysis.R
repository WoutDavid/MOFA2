#importing the seurat object and the MOFA model
seurat <- readRDS("seurat_human.RDS")
mofa <- load_model("first_model.hdf5")
#Add metadata
samples_metadata(mofa) <- seurat@meta.data %>%
  tibble::rownames_to_column("sample") %>%
  as.data.table

#Assess correlation
plot_factor_cor(mofa)

#Variance decomposition
plot_variance_explained(mofa, max_r2 = 4)
plot_variance_explained(mofa, plot_total = TRUE)[[2]]

#Characterization of factors
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
seurat <- RunPCA(seurat, npcs = K, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'pca', dims = 1:K, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes() + scale_fill_manual(values=colors)


##Using only ATAC
DefaultAssay(seurat) <- "ATAC_distal"
seurat <- RunSVD(seurat, n = K, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:K, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes() + scale_fill_manual(values=colors)












