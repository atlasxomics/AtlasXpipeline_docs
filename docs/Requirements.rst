**Requirements**
##############################
1. First, make sure you have R and RStudio installed on your computer. If you don't already have them, you can download and install R from the CRAN website (https://cran.r-project.org/) and RStudio from the RStudio website (https://rstudio.com/).

2. Once you have R and RStudio installed, you can open RStudio and create a new R script by clicking on the "File" menu and selecting "New File" followed by "R Script".

3. In your new R script, you can install ArchR by running the following command:

.. code-block::

 install.packages("ArchR")
 
4. To install Seurat, you will need to first install some other R packages that Seurat depends on. You can do this by running the following commands:
 
.. code-block::

  install.packages(c("devtools", "BiocManager"))
  BiocManager::install("Seurat")

5. Once you have installed ArchR and Seurat, you can load them into your R session by running the following commands:

.. code-block::

  library(ArchR) 
  library(Seurat)

At this point, you should be ready to use ArchR and Seurat. You can refer to the documentation for each tool to learn how to use them for your specific analysis tasks.
ArchR: https://www.archrproject.com/
Seurat: https://satijalab.org/seurat/
