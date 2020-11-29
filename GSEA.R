##################################################################################################################
## WATCH OUT: this script is reused for human/mouse seperately, so dont forget to replace all mentions of human ##
## with human or vice versa, otherwise things will go wrong.                                                    ##
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
  dcast(gs_name~gene_symbol, value.var="id", fill=as.numeric(0)) %>% 
  as.matrix %>% sapply(as.numeric)

##instead, i take the premade matrix from MOFA
#The human dataset needs a transform, because its features are written in the following form: ENSG0000000XXXXXXXX
#and our dataset has it's feature names written in regular uppercase genenames, so we use a small python script to 
#convert these identifiers, for that we write the colnames of the gene set matrix to a file
#this turned out to be harder than expected because one identifier magically dissapeared, meaning that replacing 
# the colnames in the matrix didn't fit anymore. I repeated this manually and retained some more info, which i 
#parse and then I replace the colnames.

data("MSigDB_v6.0_C5_human")

##following is commented out because it requires manual steps that I already performed, you can continue to the part where
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

#load mofa model
mofa_human <- readRDS("/media/david/Puzzles/IBP/human/second_model_human/mofa_object_human_second_model.RDS")
mofa_human

##GSEA
############## dont forget to change to the correct feature set
  
##To match the gene names,wwe have to capitalize the features
features_names(mofa_human)[["RNA"]] <- toupper(features_names(mofa_human)[["RNA"]])
features_names(mofa_human)[["RNA"]]
gsea.positive <- run_enrichment(mofa_human, 
                                feature.sets = MSigDB_v6.0_C5_human, 
                                view = "RNA",
                                sign = "positive"
)

# GSEA on negative weights
gsea.negative <- run_enrichment(mofa_human, 
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

##we see some interesting genes
genes <- list("APOA1","APOE")

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
##we see some interesting genes
genes <- list("APOA1","APOE")

genes %>% map(~ plot_factors(mofa_human, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = F
)) %>% cowplot::plot_grid(plotlist=., nrow=1)


#lets not rely only on the p-values

