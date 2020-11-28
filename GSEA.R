##################################################################################################################
## WATCH OUT: this script is reused for mouse/human seperately, so dont forget to replace all mentions of mouse ##
## with mouse or vice versa, otherwise things will go wrong.                                                    ##
##################################################################################################################
library(data.table)
library(msigdbr)
library(MOFA2)
#BiocManager::install("MOFAdata")
library(MOFAdata)

##load gene set annotations -------> the code snippet here is taken from the vignette, however it doesn't yet work because it is not 
#in the correct form for some reason, have not figure out why yet 
##[,id:=1] creates a new column with all 1's for every element
##value.var="id" chooses the values 
msgidb.matrix <- msigdbr(
  species = "Homo sapiens",
  category = "C5", 
  subcategory = "BP"
) %>% as.data.table %>% .[,id:=1] %>%
  dcast(gs_name~gene_symbol, value.var="id", fill=0) %>% 
  as.matrix

##instead, i take the premade matrix from MOFA
data("MSigDB_v6.0_C5_mouse")
head(colnames(MSigDB_v6.0_C5_mouse))

#load mofa model
mofa_mouse <- readRDS("/media/david/Puzzles/IBP/mouse/second_model_mouse/mofa_object_mouse_second_model.RDS")
mofa_mouse

##GSEA
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

##we see some interesting genes
genes <- list("APOA1","APOE")

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
##we see some interesting genes
genes <- list("APOA1","APOE")

genes %>% map(~ plot_factors(mofa_mouse, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)


#lets not rely only on the p-values

