library(Seurat)
library(dplyr)

mouse_data<- readRDS("../seurat_mouse.RDS")
mouse_data


# QC

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# mouse_data[["percent.mt"]] <- PercentageFeatureSet(mouse_data, pattern = "^MT-")
# mouse_data[["percent.mt"]] # percentage of mitochondrial genes
# 
# VlnPlot(mouse_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(mouse_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mouse_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1
# plot2
# plot1 + plot2
# 
# 
# # Keep only cells with >200 features but <2500 && < 5% mitochondrial genes
# # same numbers as the tutorial, may be worth to play witht these
# mouse_data <- subset(mouse_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Load David's preprocessed data

mouse_data <- readRDS("../third_create_seurat_mouse.RDS")


# Normalization
mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 5000)
top20 <- head(VariableFeatures(mouse_data), 20)
top20

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse_data)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#plot1 + plot2
plot1
plot2


# Scale data (mean=0, var=1 across genes in each cell)
all.genes <- rownames(mouse_data)
mouse_data <- ScaleData(mouse_data, features = all.genes)
# mouse_data[["RNA"]]@scale.data

# Dimensionality reduction

# Try PCA
mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data))

print(mouse_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouse_data, dims = 1:2, reduction = "pca")

DimPlot(mouse_data, reduction = "pca") # Bad separation
DimHeatmap(mouse_data, dims = 1, cells = 500, balanced = TRUE)


# How many PCs we should keep in PCA:
mouse_data <- JackStraw(mouse_data, num.replicate = 100)
mouse_data <- ScoreJackStraw(mouse_data, dims = 1:20)
JackStrawPlot(mouse_data, dims = 1:20)
ElbowPlot(mouse_data)

# Clustering PCA

mouse_data <- FindNeighbors(mouse_data, reduction="pca" ,dims = 1:15)
mouse_data <- FindClusters(mouse_data, resolution = 0.5)
head(Idents(mouse_data), 5)


# TSNE
mouse_data <- RunTSNE(mouse_data, dims=1:15)
DimPlot(mouse_data, reduction = "tsne", label=TRUE) # TSNE works fine enough

# Clustering TSNE

mouse_data <- FindNeighbors(mouse_data, reduction="tsne" ,dims = 1:2)
mouse_data <- FindClusters(mouse_data, resolution = 0.5)
head(Idents(mouse_data), 5)


# UMAP
#reticulate::py_install(packages ='umap-learn')

mouse_data <- RunUMAP(mouse_data, dims = 1:15)
DimPlot(mouse_data, reduction = "umap", label=TRUE) # UMAP works very well
# UMAP stored in mouse_data@reductions$umap

# Clustering UMAP

mouse_data <- FindNeighbors(mouse_data, reduction="umap" ,dims = 1:2) #k.param = 30
mouse_data <- FindClusters(mouse_data, resolution = 0.5)
head(Idents(mouse_data), 5)

# Save results so far
#saveRDS(mouse_data, file = "../seurat_mouse_preprocessed.rds")



###############################################################################
#####  Finding differentially expressed features (cluster biomarkers)   #######
###############################################################################

# find all markers of cluster 1
# cluster1.markers <- FindMarkers(mouse_data, ident.1 = 1, min.pct = 0.25) #, test.use = "roc") for DE testing
# head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(mouse_data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)


for(i in 0:21) {
  variable_name = paste("markers.cluster", i, sep = "")
  assign(variable_name, FindMarkers(mouse_data, ident.1 = i, min.pct = 0.25))
  print(paste("Head of mouse biomarkers for cluster", i, sep=" "))
  print(head(eval(as.name(variable_name)), n = 5))
}

# Results stored in markers.cluster0, markers.cluster1 ... markers.cluster21


# find markers for every cluster compared to all remaining cells, report only the positive ones
mouse_data.markers <- FindAllMarkers(mouse_data,ct = 0.25, logfc.threshold = 0.25)
mouse_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) #wt -> the variable to use for ordering results

top10 <- mouse_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mouse_data, features = top10$gene) + NoLegend()


# Re cluster data to find better clustering

mouse_data <- FindNeighbors(mouse_data, reduction="umap" ,dims = 1:2) 
mouse_data <- FindClusters(mouse_data, resolution = 0.025) 
DimPlot(mouse_data, reduction = "umap", label=TRUE) # Check clusters in UMAP plot

mouse_data.markers <- FindAllMarkers(mouse_data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
mouse_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- mouse_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mouse_data, features = top10$gene) + NoLegend()


# Write all markers
write.table(x=mouse_data.markers, file = "Markers/Mouse Markers.csv", quote = F, row.names = F, col.names = T, sep = ",")

# Write markers per cluster in ascending and descending order for Gorilla Analysis
# for (i in 0:8){
#   markers_subset = subset(mouse_data.markers, cluster==i)
#   
#   ascending = markers_subset[order(markers_subset$avg_logFC, decreasing = FALSE),]
#   write.table(x=ascending$gene, file = paste("Markers/mouse_markers_cluster", i, "_ascending.txt", sep=""), 
#               sep = ",", quote = F, col.names = F, row.names = F)
#   
#   descending = markers_subset[order(markers_subset$avg_logFC, decreasing = TRUE),]
#   write.table(x=descending$gene, file = paste("Markers/mouse_markers_cluster", i, "_descending.txt", sep=""), 
#               sep = ",", quote = F, col.names = F, row.names = F)
# }


###############################################################################
#####                           Cell annotation                         #######
###############################################################################

# https://github.com/ZJUFanLab/scCATCH
#library(devtools)
#devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)

#clu_markers <- findmarkergenes(mouse_data, species = "Mouse", cluster = 'All', match_CellMatch = FALSE, cancer = NULL,
#                                tissue = NULL, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05)
# 
# clu_ann <- scCATCH(clu_markers$clu_markers, species = "Mouse", cancer = NULL, tissue = "Brain")
# write.table(x=clu_ann,  file = "Mouse Cell Type Annotation per Cluster David.tsv", row.names = F, 
#             col.names = T, sep = "\t", quote = F)
# 
# saveRDS(clu_markers, file="scCATCH Mouse Cluster Markers.Rds")
# saveRDS(clu_ann, file = "scCATCH Mouse Cluster Annotation.Rds")

# Uncomment above lines to run scCATCH (takes a while)

clu_markers = readRDS("scCATCH Mouse Cluster Markers.Rds")
clu_ann = readRDS("scCATCH Mouse Cluster Annotation.Rds")

clu_ann$cell_type
clu_ann$celltype_score
clu_ann$cluster


# Add metadata for downstream analysis
types = clu_ann$cell_type
metadata = as.data.frame(Idents(mouse_data), row.names = rownames(mouse_data@meta.data))
colnames(metadata) = "Cluster"
# metadata$cell_type = types[metadata$Cluster + 1]

head(metadata)

metadata = metadata %>%
  mutate(cell_type = case_when(
    Cluster == 0 ~ types[1],
    Cluster == 1 ~ types[2],
    Cluster == 2 ~ types[3],
    Cluster == 3 ~ types[4],
    Cluster == 4 ~ types[5],
    Cluster == 5 ~ types[6],
    Cluster == 6 ~ types[7]
  ))

rownames(metadata) = rownames(mouse_data@meta.data)
mouse_data = AddMetaData(object = mouse_data, metadata = metadata)
head(mouse_data@meta.data$cell_type)

# Save results
saveRDS(mouse_data, file = "../Mouse_Annotated.Rds") 
