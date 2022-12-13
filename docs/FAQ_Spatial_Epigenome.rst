Frequently Asked Questions
###############################

*What is the Spatial Epigenome?*

The Spatial Epigenome is a method for analyzing the spatial organization of epigenetic modifications in a tissue sample. This can provide insights into how gene expression is regulated in different regions of a tissue, and how this regulation may be altered in disease states.


*What are the required inputs for basic downstream spatial ATAC analysis using the DBiT-seq protocol?*

The inputs needed for the Spatial Epigenome analysis include data species, number of threads, tile size, genome size, minimum TSS, minimum fragments, and set seed. These parameters can be specified when using the ArchR R package for downstream spatial ATAC analysis prepared by the DBiT-seq protocol.


*What are the dependencies for the Spatial Epigenome?*

The dependencies for the Spatial Epigenome analysis include the ArchR, Seurat, patchwork, gridExtra, kableExtra, dplyr, tibble, clusterProfiler, org.Mm.eg.db, repr, and purrr R libraries. These libraries are used for various tasks, such as creating a spatially resolved ATAC object, mapping ATAC gene scores to tissue histology, and performing functional enrichment analysis.


*How is the genome divided into smaller regions for analysis? What is the purpose of the tile size parameter?*

The tile size parameter specifies the size of the tiles that will be used to divide the genome into smaller regions for analysis. This parameter can affect the resolution of the analysis, with smaller tile sizes providing greater resolution but also requiring more computational resources.


*What is the purpose of the minimum TSS and minimum fragments parameters?*

The minimum TSS and minimum fragments parameters are used to filter out poor quality data from the analysis. These parameters specify the minimum number of transcription start sites and fragments that must be detected in a region for it to be considered significant.


*What is the purpose of the set seed parameter?*

The set seed parameter is used to specify a seed value for the random number generator used by the analysis. This can be useful for reproducing results, as it ensures that the same random number sequence will be generated each time the analysis is run with the same seed value.


*How can poor quality tixels be filtered out from the dataset in spatial ATAC analysis?*

Poor quality tixels can be filtered out from the dataset in spatial ATAC analysis by specifying minimum TSS and minimum fragment values as input parameters during the ArrowFile creation step.


*Can the ArchR ArrowFile creation step be customized to include only specific tixels in the downstream analysis?*

Yes, the ArchR ArrowFile creation step can be customized to include only specific tixels in the downstream analysis by using the subsetArchRProject function, which takes a list of tixel names as input and outputs a filtered ArchR project containing only the specified tixels.


*Can the methods and concepts demonstrated in this tutorial be applied to other spatial omics sequencing protocols or data types?*

Yes.


*What is the role of the patchwork, gridExtra, and kableExtra libraries in the Spatial Epigenome analysis?*

The patchwork, gridExtra, and kableExtra libraries are used in the Spatial Epigenome analysis to create customizable visualizations and tables of the downstream analysis results.


*What is the purpose of the clusterProfiler and org.Mm.eg.db libraries in the Spatial Epigenome analysis?*

The clusterProfiler and org.Mm.eg.db libraries are used in the Spatial Epigenome analysis to perform gene set enrichment analysis and identify marker genes for each cluster of cells.


*How does the repr and purrr libraries fit into the Spatial Epigenome analysis workflow?*

The repr and purrr libraries are used in the Spatial Epigenome analysis to improve the efficiency and simplicity of the analysis workflow by enabling the use of functional programming techniques.


*What is the significance of the D357 mouse embryo dataset in the Spatial Epigenome analysis?*

The D357 mouse embryo dataset is used in the Spatial Epigenome analysis as an example dataset to demonstrate the downstream analysis steps and visualize the resulting data.
