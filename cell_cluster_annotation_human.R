library(Seurat)
library(dplyr)

human_data<- readRDS("../seurat_human.RDS")
human_data


# QC

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# human_data[["percent.mt"]] <- PercentageFeatureSet(human_data, pattern = "^MT-")
# human_data[["percent.mt"]] # percentage of mitochondrial genes
# 
# VlnPlot(human_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(human_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(human_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1
# plot2
# plot1 + plot2
# 
# 
# # Keep only cells with >200 features but <2500 && < 5% mitochondrial genes
# # same numbers as the tutorial, may be worth to play witht these
# human_data <- subset(human_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## Load David's preprocessed data

human_data <- readRDS("../third_create_seurat_human.RDS")


# Normalization

human_data <- NormalizeData(human_data, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)

human_data <- FindVariableFeatures(human_data, selection.method = "vst", nfeatures = 5000)
top20 <- head(VariableFeatures(human_data), 20)
top20

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(human_data)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#plot1 + plot2
plot1
plot2

# Scale data (mean=0, var=1 across genes in each cell)
all.genes <- rownames(human_data)
human_data <- ScaleData(human_data, features = all.genes)
# human_data[["RNA"]]@scale.data

# Dimensionality reduction

# Try PCA
human_data <- RunPCA(human_data, features = VariableFeatures(object = human_data))

print(human_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(human_data, dims = 1:2, reduction = "pca")

DimPlot(human_data, reduction = "pca") # Quite awful separation
DimHeatmap(human_data, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(human_data, dims = 1:15, cells = 500, balanced = TRUE)

# How many PCs we should keep in PCA:
human_data <- JackStraw(human_data, num.replicate = 100)
human_data <- ScoreJackStraw(human_data, dims = 1:20)
JackStrawPlot(human_data, dims = 1:20)
ElbowPlot(human_data)

# Clustering PCA

human_data <- FindNeighbors(human_data, reduction="pca" ,dims = 1:15)
human_data <- FindClusters(human_data, resolution = 0.5)
head(Idents(human_data), 5)


# TSNE
human_data <- RunTSNE(human_data, dims=1:15)
DimPlot(human_data, reduction = "tsne", label=TRUE) # TSNE works good

# Clustering TSNE

human_data <- FindNeighbors(human_data, reduction="tsne" ,dims = 1:2)
human_data <- FindClusters(human_data, resolution = 0.5)
head(Idents(human_data), 5)


# UMAP
#reticulate::py_install(packages ='umap-learn')

human_data <- RunUMAP(human_data, dims = 1:15)
DimPlot(human_data, reduction = "umap", label=TRUE) # UMAP works very well
# UMAP stored in human_data@reductions$umap

# Clustering UMAP

human_data <- FindNeighbors(human_data, reduction="umap" ,dims = 1:2) #k.param = 30
human_data <- FindClusters(human_data, resolution = 0.5)
head(Idents(human_data), 5)

# Save results so far
#saveRDS(human_data, file = "../seurat_human_preprocessed.rds")



###############################################################################
#####  Finding differentially expressed features (cluster biomarkers)   #######
###############################################################################

# find all markers of cluster 1
cluster1.markers <- FindMarkers(human_data, ident.1 = 1, min.pct = 0.25) #, test.use = "roc") for DE testing
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(human_data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)


for(i in 0:23) {
  variable_name = paste("markers.cluster", i, sep = "")
  assign(variable_name, FindMarkers(human_data, ident.1 = i, min.pct = 0.25))
  print(paste("Head of human biomarkers for cluster", i, sep=" "))
  print(head(eval(as.name(variable_name)), n = 5))
}

# Results stored in markers.cluster0, markers.cluster1 ... markers.cluster23


# find markers for every cluster compared to all remaining cells, report only the positive ones
human_data.markers <- FindAllMarkers(human_data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
human_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) #wt -> the variable to use for ordering results

# Visualize features for cluster 4 (for example)
VlnPlot(human_data, features = c("SLC1A2", "GPC5"))

FeaturePlot(human_data, features = c("PLP1", "CRYAB", "POLR2F", "SH3TC2-DT"))#"LINC01608", "NLGN1",reduction="umap"))# 
                               #"PDE1A", "SLC5A11", "SLC1A2", "GPC5"), reduction="umap")


top10 <- human_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(human_data, features = top10$gene) + NoLegend()


# Re cluster data to find better clustering

human_data <- FindNeighbors(human_data, reduction="umap" ,dims = 1:2) 
human_data <- FindClusters(human_data, resolution = 0.0015) 
DimPlot(human_data, reduction = "umap", label=TRUE) 

human_data.markers <- FindAllMarkers(human_data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
human_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- human_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(human_data, features = top10$gene) + NoLegend()

# Write all markers
write.table(x=human_data.markers, file = "Markers/Human Markers.csv", quote = F, row.names = F, col.names = T, sep = ",")

# Write markers per cluster in ascending and descending order for Gorilla Analysis
# for (i in 0:8){
#   markers_subset = subset(human_data.markers, cluster==i)
#   
#   ascending = markers_subset[order(markers_subset$avg_logFC, decreasing = FALSE),]
#   write.table(x=ascending$gene, file = paste("Markers/human_markers_cluster", i, "_ascending.txt", sep=""), 
#               sep = ",", quote = F, col.names = F, row.names = F)
#   
#   descending = markers_subset[order(markers_subset$avg_logFC, decreasing = TRUE),]
#   write.table(x=descending$gene, file = paste("Markers/human_markers_cluster", i, "_descending.txt", sep=""), 
#               sep = ",", quote = F, col.names = F, row.names = F)
# }


###############################################################################
#####                        Cell annotation                            #######
###############################################################################

# https://github.com/ZJUFanLab/scCATCH
#library(devtools)
#devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)

# clu_markers <- findmarkergenes(human_data, species = "Human", cluster = 'All', match_CellMatch = FALSE, cancer = NULL,
#                                tissue = NULL, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05)
# 
# clu_ann <- scCATCH(clu_markers$clu_markers, species = "Human", cancer = NULL, tissue = "Brain")
# write.table(x=clu_ann,  file = "Human Cell Type Annotation per Cluster David.tsv", row.names = F, 
#             col.names = T, sep = "\t", quote = F)

#saveRDS(clu_markers, file="scCATCH Human Cluster Markers.Rds")
#saveRDS(clu_ann, file = "scCATCH Human Cluster Annotation.Rds")

#Uncomment above code to run scCATCH (it takes a while)

clu_markers = readRDS("scCATCH Human Cluster Markers.Rds")
clu_ann = readRDS("scCATCH Human Cluster Annotation.Rds")

clu_ann$cell_type
clu_ann$celltype_score
clu_ann$cluster

# Write metadata for downstream analysis
types = clu_ann$cell_type
metadata = as.data.frame(Idents(human_data), row.names = rownames(human_data@meta.data))
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
    Cluster == 6 ~ types[7],
    Cluster == 7 ~ types[8],
    Cluster == 8 ~ types[9]
  ))

rownames(metadata) = rownames(human_data@meta.data)
human_data = AddMetaData(object = human_data, metadata = metadata)
head(human_data@meta.data$cell_type)

# Save results
saveRDS(human_data, file = "../Human_Annotated.Rds")
