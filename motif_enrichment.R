##################################################################################################################
## WATCH OUT: this script is reused for human/mouse seperately, so dont forget to replace all mentions of human ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################
library(Seurat)
library(MOFA2)
library(Signac)


############################
## Motif enrichment mouse ##
############################

seurat_mouse <- readRDS("/media/david/Puzzles/IBP/mouse/second_model_mouse/seurat_mouse_second_model.RDS")
mofa_mouse <- readRDS("/media/david/Puzzles/IBP/mouse/second_model_mouse/mofa_object_mouse_second_model.RDS")

motif.matrix <- t(as.matrix(seurat_mouse[["ATAC_distal"]]@motifs@data))
rownames(motif.matrix)
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

############################
## Motif enrichment human ##
############################

seurat_human <- readRDS("/media/david/Puzzles/IBP/human/second_model_human/seurat_human_second_model.RDS")
mofa_human <- readRDS("/media/david/Puzzles/IBP/human/second_model_human/mofa_object_human_second_model.RDS")

motif.matrix <- t(as.matrix(seurat_human[["ATAC_distal"]]@motifs@data))
rownames(motif.matrix)
motif.enrichment.positive <- run_enrichment(mofa_human,
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
MotifPlot(seurat_human[["ATAC_distal"]], motifs = sig.motifs.positive)

sig.motifs.negative <- motif.enrichment.negative$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat_human[["ATAC_distal"]], motifs = sig.motifs.negative)
