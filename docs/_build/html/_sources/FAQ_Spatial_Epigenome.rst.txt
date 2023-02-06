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

*I am unable to install Seurat, what should I do?*

Make sure that you have the latest version of R installed. Seurat may not be compatible with older versions of R.
Make sure that you have the necessary dependencies installed. Seurat requires other R packages such as "dplyr", "ggplot2", "Matrix", and "SummarizedExperiment" to be installed.

*I am getting an error message about Rtools when trying to install Seurat, what should I do?*

If you are getting an error message about Rtools when trying to install Seurat, it means that Rtools is not installed on your system. Rtools is a collection of tools that are necessary for building R packages on Windows. You can download the latest version of Rtools from the CRAN website.

*I am getting an error message about Bioconductor when trying to install Seurat, what should I do?*

Seurat is available on Bioconductor, an R package repository for bioinformatics packages. If you are getting an error message about Bioconductor, it means that the Bioconductor repository is not installed on your system. You can install Bioconductor by running the following command in R: source("https://bioconductor.org/biocLite.R") then biocLite("Seurat")

*I am unable to install patchwork, what should I do?*

Make sure that you have the necessary dependencies installed. patchwork requires some other R packages such as "ggplot2" and "gridExtra" to be installed.
patchwork is available on CRAN, an R package repository. If you are getting an error message about Bioconductor, it means that you are trying to install patchwork from the wrong repository. Bioconductor is for bioinformatics packages. You should install patchwork from CRAN by running the following command in R: install.packages("patchwork")

*Why do we include certain dependencies like gridExtra, patchwork, kableExtra etc. if they only relate to plotting/graphics handling?*

These dependencies are included because they provide useful functions and features for creating and customizing plots and tables in R. The ArchR package is an R package for spatial analysis of ATAC-seq data, and creating plots and tables is an important aspect of visualizing and interpreting the results of this analysis. The dependencies like gridExtra, patchwork, kableExtra, etc. are used to create and customize plots and tables in a way that is aesthetically pleasing and easy to read. Additionally, these dependencies also provide more advanced functionality such as grid layout, table formatting and annotation, which is helpful to make the results more informative and easy to understand.

*Are these dependencies necessary for the functionality of the ArchR package?*

No, these dependencies are not necessary for the core functionality of the ArchR package. The ArchR package provides functions for performing spatial analysis of ATAC-seq data and the dependencies are not required for the analysis itself. However, they are required if you want to use the built-in functions for creating plots and tables, which is an important part of visualizing and interpreting the results of the analysis.

