# MOFA2 Analysis pipeline for ATAC and RNA-seq data
## Integrated Bioinformatics Project
### Authors: David Wouters, Tine Logghe, Jose Ignacio Alvira Larizgoitia, Alexander Ian Taylor

This github contains the code for the authors' integrated bioinformatics project, as a part of their evaluation of the eponymous course at the KU Leuven (2020-2021), as part of their 2nd year Master of Bioinformatics program. 

- **Pipeline Code**: Contains all R-scripts that were used in the analysis pipeline
  - ***create_seurat.R***: This Rscript was used to transform the CellRanger output (countmatrix + metadata) into a Seurat object. Quality control is also performed.
  - ***trainingTheModel.R***: This Rscript contains all code that adds metadata that we ourselves generated, and that was not just output of CellRanger. 
After which, the data was normalized, scaled and variable features were found. This final Seurat object is then transformed to a MOFA object. The last few lines add the MOFA model specifications and actually train the model.
  - ***downstream_analysis.R***: Contains all code that creates the plots and allows us to infer our conclusions written in the paper.
  - ***cell_cluster_annotation_organism.R***: Contains all code that was used to annotate the cells with a celltype based on UMAP clustering and their marker genes.
  - ***GSEA/motif_enrichtmen.R***: These scripts contain the code used to perform GSEA and motif enrichment on the two data modalities. When functional, this code was added to the downstream analysis, so these files can be disregarded.

 
 
 For more information, contact us at pleasedontcontactus@student.kuleuven.be
