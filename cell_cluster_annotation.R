library(Seurat)

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





