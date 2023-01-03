Spatial Transcriptome
_____________________

This program is designed to analyse and generate an analysis report in HTML format for AtlasXomics spatial datasets.

Requirements
###############

Install the following required packages and libraries using the install.packages() function. ::

  install.packages('Seurat')
  install.packages('SeuratData')
  install.packages('ggplot2')
  install.packages('patchwork')
  install.packages('dplyr')
  install.packages('devtools')
  install.packages('DropletUtils')
  install.packages('kableExtra')
  install.packages('magick')
  install.packages('stringr')
  install.packages('getopt')
  install.packages('plyr')
  
If the BiocManager package is not already installed, install it: ::

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
Use the BiocManager to install the DropletUtils and org.*.eg.db packages: ::

  BiocManager::install("DropletUtils")
  BiocManager::install("org.Hs.eg.db")
  BiocManager::install("org.Mm.eg.db")
If the BiocManager package is not already installed, install it: ::

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
Use the BiocManager to install the topGO package: ::

  BiocManager::install("topGO")
  
Load the required libraries: ::

  library(devtools)
  library(Seurat)
  library(AtlasTools)
  library(ggplot2)
  library(patchwork)
  library(plyr)
  library(DropletUtils)
  library(kableExtra)
  library(cowplot)
  library(magick)
  library(stringr)
  library(topGO)
  library(getopt)
  library(dplyr)
  library(purrr)
  library(rjson)
  library(viridis)
  library(SeuratDisk)
  
Define the command-line options and read in the values: ::

  spec = matrix (c(
    'verbose', 'v', 2, "integer",
    'help', 'h', 0, "logical",
    'folder', 'f', 1, "character",
    'species', 's', 1, "character",
    'project_name', 'n', 1, "character",
    'ref_path', 'r', 1, "character",
    'pt_size', 't', 1, "double"
  ), byrow=TRUE, ncol=4)
  
Use the getopt() function to parse the command-line options and specify the values for the necessary variables, such as the dataset folder, species, project name, reference path, and point size. ::

  opt = getopt(spec)
  dataset <- file.path(opt$folder, "Gene/raw")
  data_species <- opt$species
  project_name <- opt$project_name
  if (is.null(opt$pt_size)) {
    pt_size_factor <- 1
  } else {
    pt_size_factor <- opt$pt_size
  }
  ref_path <- opt$ref_path

Read the metadata JSON file and extract the relevant information, such as the assay type, number of barcodes, collaborator, disease state, number of channels, organ, run, species, and experiment type. Use the map() and replace() functions to handle NULL values and convert the metadata to a data frame. ::

  meta_out <- fromJSON(file = paste0(dataset, "/spatial/metadata.json")) 
  meta_rotation <- meta_out$orientation$rotation
  meta_out <- meta_out[c("assay", "barcodes", "collaborator", "diseaseState", "numChannels", "organ", "run", "species", "Experiment	type" )]
  meta_out <- map(meta_out, ~ replace(.x, is.null(.x), "NA"))
  meta_out <- as.data.frame(meta_out)

Path for flow images
########################### 

Set the path for the flow images. The fig_path variable is set to the path of the flow images, which is located in the "spatial/figure" subdirectory of the dataset directory. ::

  fig_path = file.path(dataset, "spatial/figure")
  
Check if the path for flow images exists. If it does, then iterate over the list of files to process them one by one. If the path does not exist, it creates the directory. ::

  if(dir.exists(fig_path)){
    for(file in c("flowA", "flowB", "postB_BSA", "postB\\b.", "fix")){
      skip_to_next <- FALSE
      tryCatch(
        # code to process the images goes here
      )
    }
  } else{
    dir.create(fig_path)
  }
  
Inside the loop, list the flowA.tif or flowA.TIF image. Filter it using the grep function with the file variable as the pattern. The file.path function is then used to combine the fig_path and the tif_file to get the full path of the image file. ::

  index <-  list.files(path = fig_path, pattern = "tif|TIF") %>% grep(pattern = file)
  tif_file <- list.files(path = fig_path, pattern = "tif|TIF")[index]

Read the image using the image_read function: ::

  img <- image_read(file.path(fig_path, tif_file))
  
If the file is not postB_BSA or postB\b. and the meta_rotation value is not zero, rotate the image by 360 - meta_rotation: :: 

  if((!file %in% c("postB\\b.", "postB_BSA")) && (meta_rotation !=0)){
    angle <- 360 - meta_rotation
    img <- image_rotate(img, angle)
  } 

Resize the image to a width of 950px and a height that is proportional to the width: ::

  img <- image_scale(img, "950x")

Save the image as a .png file using the image_write function. The gsub function is used to remove the extension from the tif_file variable, and the paste0 function is used to create the new filename for the .png file. ::

  base_fname <- gsub(pattern = "\\..*", replacement = "", x = tif_file)
  image_write(img, path = paste0(fig_path, "/", base_fname, "-1.png", sep=""), format = "png", quality = 75)
  
In case of any error or warning, handle them using the tryCatch function: ::

  tryCatch(
    # code to process the images goes here
    error=function(e){
      # code to handle errors goes here
    },
    warning=function(w){
      # code to handle warnings goes here
    }
  )
Create a function called Load_AtlasXomics to read the AtlasXomics spatial dataset into a Seurat object. It accepts several parameters: 

*data.dir:* The directory containing matrix.mtx, genes.tsv, barcodes.tsv along with a subdirectory spatial containing .png tissue image, scalefactors_json.json, and tissue_positions_list.csv.

*assay:* The name of the assay to assign to the data within the Seurat object (default is 'Spatial').

*slice:* The name of the tissue to assign to the data within the Seurat object (default is 'slice1').

*filter.matrix:* A logical value indicating whether to filter the spatial expression matrix based on known tissue positions (default is TRUE).

*to.upper:* A logical value indicating whether to convert the row names of the data matrix to uppercase (default is FALSE). ::

  Load_AtlasXomics <- function(
    data.dir,
    assay = 'Spatial',
    slice = 'slice1',
    filter.matrix = TRUE,
    to.upper = FALSE,
    ...
  ) {
  
Inside the function, check if the length of data.dir is greater than 1. If it is, print a warning message saying that the function only accepts one data.dir, and set data.dir equal to the first element in the list. ::

  if (length(x = data.dir) > 1) {
    warning("'Load_AtlasXomics' accepts only one 'data.dir'", immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  
Read the data using the Read10X function, passing in data.dir and ... as arguments. ::

  data <- Read10X(data.dir = data.dir, ...)
  
If to.upper is TRUE, convert the row names of data to uppercase using the toupper function. ::

  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  
Create a Seurat object from data using the CreateSeuratObject function, setting the assay and project arguments to assay and data.dir respectively. ::

  object <- CreateSeuratObject(counts = data, assay = assay, project = data.dir)
  
Set the working directory to file.path(data.dir, 'spatial') and assign the list of png files to a variable called png_file. ::

  starting_wd = getwd()
  setwd(dir = file.path(data.dir, 'spatial'))
  png_file <- list.files(pattern = "\\.png$")
  
Set the working directory back to the starting working directory. ::

  setwd(starting_wd)
  

Read the image data using the Read10X_Image function, setting the image.dir argument to file.path(data.dir, 'spatial') and the filter.matrix argument to filter.matrix. ::

  image <- Read10X_Image(
    image.dir = file.path(data.dir, 'spatial'),
    filter.matrix = filter.matrix
  )
  
Subset the image data using the Cells function and the Seurat object, and set the default assay of the image data to assay. ::

  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
Add the image data to the Seurat object using the [[]] operator, setting the slice to image.::

  object[[slice]] <- image

If filter.matrix is TRUE, subset the Seurat object to only include rows that are on tissue cells using the rownames function. ::

  if (filter.matrix) {
    on_tissue_cells <- rownames(object@images$slice1@coordinates)
    object <- object[,on_tissue_cells]
  }
  return(object)
  
Create a Seurat object called object_AXOSpatial_seurat by calling the Load_AtlasXomics function with the data.dir argument set to dataset. ::

  object_AXOSpatial_seurat = Load_AtlasXomics(data.dir = dataset)
  
Create a Seurat object called object_AXOSpatial_seurat_all_tixels by calling the Load_AtlasXomics function with the data.dir argument set to dataset and the filter.matrix argument set to FALSE. ::

  object_AXOSpatial_seurat_all_tixels = Load_AtlasXomics(data.dir = dataset, filter.matrix = FALSE)

Add a new metadata column called orig.ident to both object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels, setting the value to project_name as a factor. ::

  object_AXOSpatial_seurat$orig.ident = as.factor(project_name)
  object_AXOSpatial_seurat_all_tixels$orig.ident = as.factor(project_name)

Set the Idents of both object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels to 'orig.ident'. ::

  Idents(object_AXOSpatial_seurat) = 'orig.ident'
  Idents(object_AXOSpatial_seurat_all_tixels) = 'orig.ident'

Create a new object called object_AXOSpatial_seurat_all_tixels0 that is a copy of object_AXOSpatial_seurat_all_tixels. ::
  
  object_AXOSpatial_seurat_all_tixels0 <- object_AXOSpatial_seurat_all_tixels

Add new metadata columns to object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels using the PercentageFeatureSet function. ::

  object_AXOSpatial_seurat[["percent.mt"]] <- PercentageFeatureSet(object_AXOSpatial_seurat, pattern = "^[Mm][Tt]-")
  object_AXOSpatial_seurat[["percent.rb"]] <- PercentageFeatureSet(object_AXOSpatial_seurat, pattern = "^[Rr][Pp][Ss]|[Rr][Pp][Ll]")
  object_AXOSpatial_seurat[["percent.hg"]] <- PercentageFeatureSet(object_AXOSpatial_seurat, "^[Hh][Bb][^Pp]")

Calculate the percentage of features in the object_AXOSpatial_seurat_all_tixels object that match the pattern "^[Mm][Tt]-". Assign the result to a new slot called "percent.mt" in the object_AXOSpatial_seurat_all_tixels object. ::

  object_AXOSpatial_seurat_all_tixels[["percent.mt"]] <- PercentageFeatureSet(object_AXOSpatial_seurat_all_tixels, pattern = "^[Mm][Tt]-")
Calculate the percentage of features in the object_AXOSpatial_seurat_all_tixels object that match the pattern "^[Rr][Pp][Ss]|[Rr][Pp][Ll]". Assign the result to a new slot called "percent.rb" in the object_AXOSpatial_seurat_all_tixels object. ::

    object_AXOSpatial_seurat_all_tixels[["percent.rb"]] <-      PercentageFeatureSet(object_AXOSpatial_seurat_all_tixels, pattern = "^[Rr][Pp][Ss]|[Rr][Pp][Ll]")

Calculate the percentage of features in the object_AXOSpatial_seurat_all_tixels object that match the pattern "^[Hh][Bb][^Pp]". Assign the result to a new slot called "percent.hg" in the object_AXOSpatial_seurat_all_tixels object. ::

  object_AXOSpatial_seurat_all_tixels[["percent.hg"]] <- PercentageFeatureSet(object_AXOSpatial_seurat_all_tixels, "^[Hh][Bb][^Pp]")
Add metadata from the slice1 slot of the images slot of the object_AXOSpatial_seurat_all_tixels object to the object_AXOSpatial_seurat_all_tixels object. ::

  object_AXOSpatial_seurat_all_tixels = AddMetaData(object_AXOSpatial_seurat_all_tixels,  object_AXOSpatial_seurat_all_tixels@images$slice1@coordinates)
Add metadata from the slice1 slot of the images slot of the object_AXOSpatial_seurat object to the object_AXOSpatial_seurat object. ::

  object_AXOSpatial_seurat = AddMetaData(object_AXOSpatial_seurat, object_AXOSpatial_seurat@images$slice1@coordinates)
Create a file name based on the project_name variable and the current date and time. Save the object_AXOSpatial_seurat_all_tixels object to a file with this file name and the extension rds. ::

  obj_name = paste0(sub('\\..*', '', project_name), format(Sys.time(),'_%Y%m%d_%H%M%S'))
  saveRDS(object_AXOSpatial_seurat_all_tixels, paste0(obj_name, "_unfiltered_all_tixel.rds"))

Define a function qc_stats_df that takes in two arguments: seurat_object and row_name.
Inside the function, create a data frame df with columns: Num_Tixels, UMI_Average, UMI_std, UMI_min, UMI_max, Genes_Average, Genes_std, Genes_min, and Genes_max.
For each column, calculate the appropriate statistic (e.g. sum, mean, standard deviation, minimum, maximum) from the seurat_object and assign it to the corresponding column in df.
Set the row.names of df to row_name.
Transpose df using the t() function and return the result. ::

  qc_stats_df <- function(seurat_object, row.name){
    df <- data.frame(
      Num_Tixels = sum(seurat_object@images$slice1@coordinates$tissue),
      UMI_Average = mean(seurat_object@meta.data$nCount_Spatial),
      UMI_std = sd(seurat_object@meta.data$nCount_Spatial),
      UMI_min = min(seurat_object@meta.data$nCount_Spatial),
      UMI_max = max(seurat_object@meta.data$nCount_Spatial),
      Genes_Average = mean(seurat_object@meta.data$nFeature_Spatial),
      Genes_std = sd(seurat_object@meta.data$nFeature_Spatial),
      Genes_min = min(seurat_object@meta.data$nFeature_Spatial),
      Genes_max = max(seurat_object@meta.data$nFeature_Spatial),
      row.names = row.name
    ) %>% t()
    return(df)
  }
  
Compute the summary statistics for the object_AXOSpatial_seurat object and assign the resulting data frame to on_tiss_before_filter.
::
  on_tiss_before_filter <- qc_stats_df(object_AXOSpatial_seurat, row.name = 'Before Filtering')

Extract the meta.data from object_AXOSpatial_seurat and add 1 to the col and row columns. Assign the resulting data frame to meta.data_plots..::

  meta.data_plots = object_AXOSpatial_seurat@meta.data
  meta.data_plots[,c('col', 'row')] = meta.data_plots[,c('col', 'row')] + 1

Using the ggplot2 package, create a bar plot of UMI counts by column using meta.data_plots. Label the y-axis as "#UMIs / row" and the x-axis as "Row (project_name)". Store the plot in a variable umi_row. ::


  umi_row <-ggplot(data=meta.data_plots, aes(x=col, y=nCount_Spatial)) + geom_bar(stat="identity") + ylab('#UMIs / row') + xlab(paste('Row (', project_name, ')', sep="")) + theme(text=element_text(size=21))

Create a bar plot of the nCount_Spatial values by row using ggplot2. Assign the resulting plot to umi_column.::

  umi_column <-ggplot(data=meta.data_plots, aes(x=row, y=nCount_Spatial)) +
    geom_bar(stat="identity") + ylab('#UMIs / column') + xlab(paste('Column (', project_name, ')', sep="")) + theme(text=element_text(size=21))

  umi_QC = wrap_plots(umi_row, umi_column)

Find the 5 rows and 5 columns with the lowest UMI counts.

Aggregate the meta.data_plots$nCount_Spatial column by the meta.data_plots$row column and sum the values.
Sort the resulting data frame by the aggregated values in ascending order.
Select the first 5 rows of the sorted data frame and store it in lowest_col_UMI.
Add a new column called UMI count to lowest_col_UMI with the aggregated values.
Remove the original aggregated values column from lowest_col_UMI.
Remove row names from lowest_col_UMI.
Aggregate the meta.data_plots$nCount_Spatial column by the meta.data_plots$col column and sum the values.
Sort the resulting data frame by the aggregated values in ascending order.
Select the first 5 rows of the sorted data frame and store it in lowest_row_UMI.
Add a new column called UMI count to lowest_row_UMI with the aggregated values.
Remove the original aggregated values column from lowest_row_UMI.
Remove row names from lowest_row_UMI.
Write lowest_col_UMI to a CSV file at the file path dataset/spatial/lowest_col_UMI.csv.
Write lowest_row_UMI to a CSV file at the file path dataset/spatial/lowest_row_UMI.csv. ::

  aggregate_col <- aggregate(x = meta.data_plots$nCount_Spatial, by = list(Col = meta.data_plots$row), FUN=sum)
  aggregate_col <- aggregate_col[order(aggregate_col$x),]
  lowest_col_UMI <- aggregate_col[1:5,]
  lowest_col_UMI$`UMI count` <- lowest_col_UMI$x
  lowest_col_UMI$x <- NULL
  rownames(lowest_col_UMI) <- NULL
  aggregate_row <- aggregate(x = meta.data_plots$nCount_Spatial, by = list(Row = meta.data_plots$col), FUN=sum)
  aggregate_row <- aggregate_row[order(aggregate_row$x),]
  lowest_row_UMI <- aggregate_row[1:5,]
  lowest_row_UMI$`UMI count` <- lowest_row_UMI$x
  lowest_row_UMI$x <- NULL
  rownames(lowest_row_UMI) <- NULL
  write.csv(x=lowest_col_UMI, file=file.path(dataset, 'spatial', 'lowest_col_UMI.csv'))
  write.csv(x=lowest_row_UMI, file=file.path(dataset, 'spatial', 'lowest_row_UMI.csv'))

Filter object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels to only include rows where nFeature_Spatial is greater than 50 and percent.mt is less than 30. ::

  object_AXOSpatial_seurat <- subset(object_AXOSpatial_seurat, subset = nFeature_Spatial > 50 & percent.mt < 30)

For object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels:

Convert the Spatial@counts assay data to a data frame.
Select only the columns in the data frame that have at least one non-zero value.
Filter the original object_AXOSpatial_seurat or object_AXOSpatial_seurat_all_tixels object to only include the cells that are present in the filtered data frame. ::

  object_AXOSpatial_seurat_all_tixels <- subset(object_AXOSpatial_seurat_all_tixels, subset = nFeature_Spatial > 50 & percent.mt < 30)

  mat = as.data.frame(object_AXOSpatial_seurat@assays$Spatial@counts)

  mat2 = mat[, colSums(mat != 0) > 0]

  object_AXOSpatial_seurat = subset(object_AXOSpatial_seurat, cells = colnames(mat2))

  mat = as.data.frame(object_AXOSpatial_seurat_all_tixels@assays$Spatial@counts)

  mat2 = mat[, colSums(mat != 0) > 0]

  object_AXOSpatial_seurat_all_tixels = subset(object_AXOSpatial_seurat_all_tixels, cells = colnames(mat2))

Iterate over the elements in the list c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.rb", "percent.hg"). For each element i, do the following:

Create a plot using the VlnPlot function, with object_AXOSpatial_seurat as the input and i as the features argument.
Add a box plot to the plot using the geom_boxplot function.
Remove the legend from the plot using the NoLegend function.
Assign the resulting plot to the element in a with the name i. ::

  a = c()
  #  VlnPlot(object_AXOSpatial_seurat, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.rb", "percent.hg"), pt.size = 0.0, combine = FALSE) + NoLegend()
  for(i in c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.rb", "percent.hg")){
    a[[i]] <-VlnPlot(object_AXOSpatial_seurat, features = i, pt.size = 0.0) + geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) + NoLegend()
  }
  
Combine the plots in a into a single plot using the CombinePlots function, with a layout of 5 columns.::

  a = CombinePlots(a, ncol = 5)

Create a scatter plot of nCount_Spatial vs nFeature_Spatial for object_AXOSpatial_seurat using the FeatureScatter function.

Set the point size to 1 and the color to black.
Remove the legend from the plot using the NoLegend function.
Set the text size to 21 using the theme function.
Store the resulting plot in a variable b. ::

  b = FeatureScatter(object_AXOSpatial_seurat, "nCount_Spatial", "nFeature_Spatial", pt.size = 1, cols = 'black', ) + NoLegend() +theme(text=element_text(size=21))
  
Create an empty plot using the ggdraw function and store it in a variable qcPlot. Add the plots a and b to the qcPlot plot in the specified positions and sizes using the draw_plot function. ::

  qcPlot = ggdraw() +
    draw_plot(a, x = 0, y = 1/2, width = 1, height = 1/2) +
    draw_plot(b, x = 0, y = 0, width = 1, height = 1/2)
    
Set the scientific notation threshold to 999 using the options function. Calculate statistics for object_AXOSpatial_seurat using the qc_stats_df function and store the resulting data frame in a variable on_tiss_after_filter. Bind the data frame on_tiss_before_filter and on_tiss_after_filter together horizontally and store the result in a variable on_tiss_stats.::

  options(scipen = 999)
  on_tiss_after_filter <- qc_stats_df(object_AXOSpatial_seurat, row.name = "After Filtering")
  on_tiss_stats <- cbind(on_tiss_before_filter , on_tiss_after_filter)

Set the Idents of object_AXOSpatial_seurat_all_tixels to 'tissue'. ::

  Idents(object_AXOSpatial_seurat_all_tixels) = 'tissue'

Check if the value 0 is not present in the Idents of object_AXOSpatial_seurat_all_tixels. If it is not present, do the following:
Set the Idents of object_AXOSpatial_seurat_all_tixels to 'orig.ident'.
Create a data frame with 0 values for UMI_Average_NT, UMI_std_NT, UMI_min_NT, UMI_max_NT, and percent_umi_off_tissue and a row name of 'After Filtering'.
Transpose the data frame and store it in a variable off_tiss_filtered_df. ::

  f(!0 %in% Idents(object_AXOSpatial_seurat_all_tixels)){
    Idents(object_AXOSpatial_seurat_all_tixels) = 'orig.ident'
    off_tiss_filtered_df <- data.frame(
      UMI_Average_NT = 0,
      UMI_std_NT = 0,
      UMI_min_NT = 0,
      UMI_max_NT = 0,
      percent_umi_off_tissue = 0,
      row.names = 'After Filtering'
    ) %>% t()
  } else {

If the value 0 is present in the Idents of object_AXOSpatial_seurat_all_tixels, do the following:
Filter object_AXOSpatial_seurat_all_tixels to only include cells with an ident of 0 and store the result in object_AXOSpatial_seurat_Non_Tissue.
Set the Idents of object_AXOSpatial_seurat_all_tixels to 'orig.ident'.
Calculate the mean, standard deviation, minimum, and maximum of the nCount_Spatial values for object_AXOSpatial_seurat_Non_Tissue and the percentage of UMI counts that are from off-tissue cells.
Create a data frame with these calculated values and a row name of 'After Filtering'.
Transpose the data frame and store it in a variable off_tiss_filtered_df. ::

  object_AXOSpatial_seurat_Non_Tissue = subset(object_AXOSpatial_seurat_all_tixels, ident = 0)
    Idents(object_AXOSpatial_seurat_all_tixels) = 'orig.ident'
    off_tiss_filtered_df <- data.frame(
      UMI_Average_NT = mean(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial),
      UMI_std_NT = sd(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial),
      UMI_min_NT = min(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial),
      UMI_max_NT = max(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial),
      percent_umi_off_tissue = sum(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial)/(sum(object_AXOSpatial_seurat_Non_Tissue@meta.data$nCount_Spatial)+sum(object_AXOSpatial_seurat@meta.data$nCount_Spatial)),
      row.names = 'After Filtering'
    ) %>% t()
  }

Read a CSV file at the file path dataset/../Summary.csv and store the result in a data frame Sout. ::

  Sout <- read.csv(file.path(dataset, "../Summary.csv"), header=FALSE, sep=",")

Create a data frame with the following columns and values:

run: project_name
Num_tixels: the rounded value of on_tiss_after_filter['Num_Tixels']
UMI_perc_off_tixel: the rounded value of off_tiss_filtered_df['percent_umi_off_tissue'] multiplied by 100
Gene_per_tixel: the rounded value of on_tiss_after_filter['Genes_Average']
UMI_per_tixel: the rounded value of on_tiss_after_filter['UMI_Average']
Number_of_reads: the value in the second column of Sout ::

  json_data_frame <- data.frame(run=project_name,
                                Num_tixels=round(on_tiss_after_filter['Num_Tixels'], 2),
                                UMI_perc_off_tixel=round(off_tiss_filtered_df['percent_umi_off_tissue'], 2)*100,
                                Gene_per_tixel=round(on_tiss_after_filter['Genes_Average'], 2),
                                UMI_per_tixel=round(on_tiss_after_filter['UMI_Average'], 2),
                                Number_of_reads=Sout$V2[1])

Convert the data frame to a JSON object using the toJSON function and set the indentation level to 4. Write the JSON object to a file at the file path dataset/spatial/stats.json. ::

  jsonData <- toJSON(json_data_frame, indent=4)

  write(jsonData, file.path(dataset, 'spatial', 'stats.json'))
  
Report 
#################

Calculate the number of valid reads by multiplying the value in the second column of Sout by the value in the third column of Sout.::

  n_valid_reads <- Sout$V2[1] * Sout$V2[2]
  
Add a column called log_10_nUMI to object_AXOSpatial_seurat and object_AXOSpatial_seurat_all_tixels by taking the log base 10 of the nCount_Spatial values. ::

  object_AXOSpatial_seurat$log_10_nUMI = log10(object_AXOSpatial_seurat@meta.data$nCount_Spatial)
  object_AXOSpatial_seurat_all_tixels$log_10_nUMI = log10(object_AXOSpatial_seurat_all_tixels@meta.data$nCount_Spatial) 
  
Create two plots of spatial data:

*plot_ALL:* a plot of all tixels in object_AXOSpatial_seurat_all_tixels with the log_10_nUMI values as the color, using the SpatialFeaturePlot function. The plot should have a title of "All Tixels (project_name)" and the legend position should be set to "right".

*plot_Tissue:* a plot of only the tixels on tissue in object_AXOSpatial_seurat with the log_10_nUMI values as the color, using the SpatialFeaturePlot function. The plot should have a title of "On Tissue Tixels (project_name)" and the legend position should be set to "right". 

Combine the two plots plot_ALL and plot_Tissue into a single plot using the wrap_plots function and store the result in a variable before_after_umiPlot. ::

  plot_ALL <- SpatialFeaturePlot(object_AXOSpatial_seurat_all_tixels, features = "log_10_nUMI", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right") + ggtitle(paste("All Tixels (", project_name, ")", sep="")) + labs(color = 'log10_nUMI') + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16))
  
  plot_Tissue <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "log_10_nUMI", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right") + ggtitle(paste("On Tissue Tixels (", project_name, ")", sep="")) + labs(color = 'log10_nUMI') + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16))

  before_after_umiPlot = wrap_plots(plot_ALL, plot_Tissue)

Create two plots:

plot1: a boxplot of the log_10_nUMI values for object_AXOSpatial_seurat, using the VlnPlot function.
plot2: a plot of the spatial data for object_AXOSpatial_seurat with the log_10_nUMI values as the color, using the SpatialFeaturePlot function. ::

   plot1 <- VlnPlot(object_AXOSpatial_seurat, features = "log_10_nUMI", pt.size = 0) + geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) + NoLegend() + xlab("") + ylab("#UMIs / pixel") + labs(color = 'log10_nUMI')
  plot2 <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "log_10_nUMI", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right") + labs(color = 'log10_nUMI')
  plot2

  umiPlot = wrap_plots(plot1, plot2)

Create plots that displays the number of genes per pixel for the object_AXOSpatial_seurat object. Plot1.5 is a violin plot showing the distribution of the number of genes per pixel. This plot is created using the VlnPlot() function from the seurat package and the geom_boxplot() function from the ggplot2 package. The VlnPlot() function plots the distribution of a feature for a given seurat object, and the geom_boxplot() function adds a boxplot to the plot created by VlnPlot().

plot2.5 is a spatial feature plot showing the number of genes per pixel across all pixels. This plot is created using the SpatialFeaturePlot() function from the seurat package. This function plots the distribution of a feature for a given seurat object on a 2D grid, with each grid cell corresponding to a pixel. ::

  plot1.5 <- VlnPlot(object_AXOSpatial_seurat, features = "nFeature_Spatial", pt.size = 0)  + geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) + NoLegend() + xlab("") + ylab("#Genes / pixel")
  plot2.5 <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "nFeature_Spatial", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right")

  genePlot = wrap_plots(plot1.5, plot2.5)

Create plots that displays the percentage of mitochondrial and ribosomal genes per pixel for the object_AXOSpatial_seurat object. plot3, a boxplot, shows the distribution of the percentage of mitochondria per pixel. plot4, a scatterplot, shows the percentage of mitochondria per pixel for each pixel. These plots are created using the same process as the gene plot. ::

  plot3 <- VlnPlot(object_AXOSpatial_seurat, features = "percent.mt", pt.size = 0)  + geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) + NoLegend() + xlab("") + ylab("mt% / pixel")
  plot4 <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "percent.mt", alpha = c(0.8, 1), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right")
  mtPlot = wrap_plots(plot3, plot4)

plot5, a boxplot, shows the distribution of the percentage of ribosomal genes per pixel. This plot is created using the same process above. plot6, a scatterplot, shows the percentage of ribosomal genes per pixel for each pixel. ::

  plot5 <- VlnPlot(object_AXOSpatial_seurat, features = "percent.rb", pt.size = 0)  + geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) + NoLegend() + xlab("") + ylab("rb% / pixel")
  plot6 <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "percent.rb", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) + ggplot2::theme(legend.position = "right")
  rbPlot = wrap_plots(plot5, plot6)
  
Create Plot 7 by creating a violin plot of the percent.hg feature using the VlnPlot function and set the pt.size parameter to 0. Use the geom_boxplot function to add a boxplot to the plot, with a width of 0.1, a black outline, and a white fill. Use the NoLegend function to remove the legend and the xlab and ylab functions to label the x-axis and y-axis, respectively. ::

  plot7 <- VlnPlot(object_AXOSpatial_seurat, features = "percent.hg", pt.size = 0)  +
           geom_boxplot(width=0.1, color="black", fill="white", outlier.shape = NA) +
           NoLegend() + xlab("") + ylab("hg% / pixel")
           
Create a spatial feature plot of the percent.hg feature using the SpatialFeaturePlot function and set the alpha parameter to a range of values and the pt.size.factor parameter to the pt_size_factor variable. Use the ggplot2::theme function to set the legend position to "right". ::

  plot8 <- SpatialFeaturePlot(object_AXOSpatial_seurat, features = "percent.hg", alpha = c(0.8, 2), pt.size.factor = pt_size_factor) +
           ggplot2::theme(legend.position = "right")
           
Use the wrap_plots function to combine the violin plot and spatial feature plot into a single plot called hgPlot.::

  hgPlot = wrap_plots(plot7, plot8)
  
Use the GetAssayData function to get the assay data for the Spatial assay and store it in the counts variable.::

  counts <- GetAssayData(object_AXOSpatial_seurat, assay = "Spatial")

Use the grep function to find the rows in counts that match the pattern "^[Mm][Tt]-" and store their indices in the mt.index variable. ::

  mt.index <- grep(pattern = "^[Mm][Tt]-", x = rownames(counts), value = FALSE)

Remove the rows corresponding to the mt.index from counts and update the object_AXOSpatial_seurat object with the updated counts matrix using the subset function.::

  counts <- counts[-mt.index, ]
  object_AXOSpatial_seurat <- subset(object_AXOSpatial_seurat, features = rownames(counts))

Set the Idents of the object_AXOSpatial_seurat object to 'tissue'.

  Idents(object_AXOSpatial_seurat) = 'tissue'
