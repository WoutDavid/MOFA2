#Dont judge me for this
#https://drive.google.com/drive/folders/1cXiEQihecVUAsT2Op63_MUkJN9sAZLke
#https://github.com/bioFAM/MOFA
#https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
#https://satijalab.org/seurat/vignettes.html
#https://github.com/WoutDavid/MOFA2/network/dependencies
#https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html


library(Seurat)
library(dplyr)

#seurat_mofa_human <- readRDS("../seurat_human.RDS")
pbmc<- readRDS("../seurat_human.RDS")
pbmc


# QC

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.mt"]] # percentage of mitochondrial genes

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
plot1 + plot2


# Keep only cells with >200 features but <2500 && < 5% mitochondrial genes
# same numbers as the tutorial, may be worth to play witht these
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization

# method "LogNormalize" that normalizes the feature expression measurements for each cell by 
# the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms 
# the result. Normalized values are stored in pbmc[["RNA"]]@data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(pbmc), 20)
top20

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#plot1 + plot2
plot1
plot2

# Scale data (mean=0, var=1 across genes in each cell)
# results in: pbmc[["RNA"]]@scale.data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# pbmc[["RNA"]]@scale.data

# Dimensionality reduction

# Try PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") # Quite awful clustering
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# How many PCs we should keep in PCA:
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

# Clustering PCA

pbmc <- FindNeighbors(pbmc, reduction="pca" ,dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

# UMAP
#reticulate::py_install(packages ='umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:15) # With 1:10 works also well
DimPlot(pbmc, reduction = "umap", label=TRUE) # UMAP works very very well
# UMAP stored in pbmc@reductions$umap

# Clustering UMAP

pbmc <- FindNeighbors(pbmc, reduction="umap" ,dims = 1:2) #k.param = 30
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

# TSNE
pbmc <- RunTSNE(pbmc, dims=1:15)
DimPlot(pbmc, reduction = "tsne", label=TRUE) # TSNE works also amazingly well

# Clustering TSNE

pbmc <- FindNeighbors(pbmc, reduction="tsne" ,dims = 1:2)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

# Save results so far
#saveRDS(pbmc, file = "../seurat_human_preprocessed.rds")



###############################################################################
#####  Finding differentially expressed features (cluster biomarkers)   #######
###############################################################################

## THERE ARE 15 CLUSTERS IN UMAP KNN GRAPH WITH DEFAULT K

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25) #, test.use = "roc") for DE testing
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)


for(i in 0:14) {
  variable_name = paste("markers.cluster", i, sep = "")
  assign(variable_name, FindMarkers(pbmc, ident.1 = i, min.pct = 0.25))
  print(paste("Head of biomarkers for cluster", i, sep=" "))
  print(head(eval(as.name(variable_name)), n = 5))
}

# Results stored in markers.cluster0, markers.cluster1 ... markers.cluster14


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) #wt -> the variable to use for ordering results

# Visualize features for cluster 4 (for example)
VlnPlot(pbmc, features = c("SLC1A2", "GPC5"))

FeaturePlot(pbmc, features = c("PLP1", "CRYAB", "POLR2F", "SH3TC2-DT"))#"LINC01608", "NLGN1",reduction="umap"))# 
                               #"PDE1A", "SLC5A11", "SLC1A2", "GPC5"), reduction="umap")
# Uncomment above to generate the plot with more features (but it is very big)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# We found 15 clusters, from 0 to 14, in the data. All the marker genes are stored under
# the markers.clusterX variable, where X is the cluster number. But for simplicity we list here the top3
# marker genes found in each cluster:
#   print(pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC), n=45)
#   Cluster 0: DBNDD2     PLP1      CRYAB
#   Cluster 1: POLR2F     SEC14L5   SH3TC2-DT 
#   Cluster 2: LINC01608  NLGN1     AC012494.1
#   Cluster 3: PDE1A      SLC5A11   LINC00609    
#   Cluster 4: LINC00499  SLC1A2    GPC5      
#   Cluster 5: PLP1       MAN2A1    LANCL1    
#   Cluster 6: DBNDD2     CRYAB     PLP1      
#   Cluster 7: MT-ND4     MT-CO2    MT-ND1    
#   Cluster 8: CD83       ST6GAL1   LRMDA     
#   Cluster 9: ATP1B3     DPP10     AKAP12    
#   Cluster 10: CST3      SLC1A2    GPC5      
#   Cluster 11: ROBO2     GRIK1     RBFOX1    
#   Cluster 12: TFRC      CLU       PDK4      
#   Cluster 13: GALNTL6   CNTN5     CSMD1     
#   Cluster 14: LUZP2     PCDH15    LHFPL3    
  

