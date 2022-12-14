Spatial Plots of Enriched Motifs
__________________________________

We use the ArchR package to perform motif enrichment analysis on the dataset and identify enriched motifs. We then use the Seurat package to add spatial data to the analysis and plot the spatial data using the enriched motifs as the features. This allows us to visualize which motifs are enriched in specific regions of the tissue, visualize it's spatial distribution, and gain insights into the regulation of gene expression in the tissue.


**Call peaks and add group coverages and reproducible peak sets**
##############################################################################

Use the addGroupCoverages function to call on the input ArchR project object, proj_in_tissue, to add group coverages to the object. The groupBy parameter specifies which metadata column to group the coverages by, in this case 'Clusters'. 
The addGroupCoverages function calculates the average coverage of each peak in a single-cell RNA-seq dataset, grouped by 'Clusters'. The resulting object contains the average coverage of each peak in each group, along with metadata about the peaks and the groups.::


  out_list <- tryCatch(expr = {
    proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters")


Use the findMacs2 function to find the path to the macs2 program on the system. Then, use the addReproduciblePeakSet function to call on the proj_in_tissue to add reproducible peak sets to the object. The groupBy, pathToMacs2, and genomeSize parameters are used to specify the metadata column to group the peaks by, the path to the macs2 program, and the size of the genome. The force parameter is set to TRUE to force re-running the peak calling even if it has already been performed. ::

    pathToMacs2 <- findMacs2()
    proj_in_tissue <- addReproduciblePeakSet(
      ArchRProj = proj_in_tissue,
      groupBy = "Clusters",
      pathToMacs2 = pathToMacs2,
      genomeSize = genomeSize,
      force = TRUE
   )

Add peak matrices
#######################################

Use the addPeakMatrix function to call on proj_in_tissue to add a peak matrix to the object. This matrix is used to store the peak calls.::

  proj_in_tissue <- addPeakMatrix(proj_in_tissue)


Add deviation matrices, motif enrichment
#########################################

Check if the Motif column is not in the names of the peak annotation data frame. If this is the case, then add motif annotations to the object using the addMotifAnnotations function. The motifSet parameter is set to cisbp for human and mouse datasets, and encode for all other species. The name parameter is set to Motif, and the force parameter is set to TRUE to force re-running the motif enrichment analysis even if it has already been performed. ::

  if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
    if (data_species == "hg38" || data_species == "mm10") {
      proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
    } else {
      proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "encode", name = "Motif", force = TRUE, species = getGenome(ArchRProj = proj_in_tissue))
    }
  }
 
Use the addBgdPeaks() function to add background peak information. This function takes the ArchRProj object as an input, along with the force argument, which is set to TRUE so that it'll overwrite any existing background peak information in the object. ::

  proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)
 
Add a matrix of deviations using the addDeviationsMatrix() function. This function takes the ArchRProj object as an input, along with the peakAnnotation argument, which specifies the name of the peak annotations to use when calculating the deviations. ::
 
  proj_in_tissue <- addDeviationsMatrix(
      ArchRProj = proj_in_tissue, 
      peakAnnotation = "Motif",
      force = TRUE
    )
  
Save project as RDS file
#######################################
Save the project as an RDS file using the saveRDS() function. RDS files are a binary file format so it can be loaded and used in future analyses ::

  saveRDS(proj_in_tissue, paste0(project_name, "_spatial_markerMotifs.rds"))

Get marker features and create list of enriched motifs
##############################################################################

Use getMarkerFeatures() to identify marker features within the ArchRProj object. The identified markers are then filtered using getMarkers() and stored in the motifs variable. ::

  markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
  )

  markerMotifsList <- getMarkers(markersMotifs,
  motifs <- list()
    for (i in seq_len(length(markerMotifsList))) {
      if (length(markerMotifsList[[i]]$name)>1) {
        motifs <- c(motifs, markerMotifsList[[i]]$name[[1]])
        motifs <- c(motifs, markerMotifsList[[i]]$name[[2]])
      }
    }
      
   
If the input list of motifs has more than one element, converts the motif to a string, and add a "z:" prefix to each motif, remove duplicate motifs, and assign the resulting list of motifs to the variable motifs. ::

     if (length(motifs)>1) {
       motifs <- unlist(motifs)
       motifs <- paste0('z:', motifs)
   motifs <- unique(motifs)


Apply addImputeWeights to the input Seurat object and assign the result to the variable proj_in_tissue. ::

  proj_in_tissue <- addImputeWeights(proj_in_tissue)

Deviation scores and matrices
#####################################

Apply getDeviation_ArchR to the modified Seurat object and the list of motifs, along with the result of applying the getImputeWeights function to the modified Seurat object. Assign the result to the variable dev_score. ::

  dev_score <- getDeviation_ArchR(ArchRProj = proj_in_tissue, name = motifs, imputeWeights = getImputeWeights(proj_in_tissue))

Set all NA values in dev_score to 0. ::

  dev_score[is.na(dev_score)] <- 0 #min(dev_score, na.rm = TRUE)

Create a new Seurat object using the dev_score matrix and the metadata from the input Seurat object, and assign the result to the variable object. ::

  object <- CreateSeuratObject(counts = dev_score, assay = "Spatial", meta.data = meta.data)
  
Filtering and setting default assay
######################################

Load image from a specified directory, filter the image based on the cells present in the object Seurat object, and set the image as the default assay for object.

Assign object to the variable spatial.obj. ::

  image <- Read10X_Image(image.dir = spatialFolder, filter.matrix = TRUE)
      image <- image[Cells(x = object)]
      DefaultAssay(object = image) <- "Spatial"
      object[['slice1']] <- image

  spatial.obj <- object

Creating Spatial plots for enriched motifs
################################################

Create a list of plots called motif_list. For each enriched motif in the spatial.obj object, create a plot using SpatialPlot_new(). The features argument specifies the motif to plot, and the pt.size.factor argument specifies the size of the points on the plot. The image.alpha and stroke arguments control the transparency and stroke width of the plot. The alpha argument controls the transparency of the points on the plot. The min.cutoff and max.cutoff arguments specify the minimum and maximum values to include on the plot. Then sets the shape of the points to squares using the shape parameter. Add the resulting plot to motif_list. ::

  motif_list <- list()
      for(i in rownames(x=spatial.obj)){
        motif_list[[i]] <- SpatialPlot_new(spatial.obj, features=i, pt.size.factor = pt_size_factor, 
                                           image.alpha = 0, stroke = 0, alpha = c(1, 1),  min.cutoff = "q10", max.cutoff = "q90") + 
          theme(legend.position = "top", legend.text=element_text(size=9), legend.title=element_text(size=9))
        motif_list[[i]]$layers[[1]]$aes_params <- c(motif_list[[i]]$layers[[1]]$aes_params, shape=22) # set spots to square shape 
      }
  
Create a combined plot of all the individual motif plots using the wrap_plots function, specifying the number of columns. ::

  motif_plots <- wrap_plots(motif_list, ncol = 3)
    
Save the combined plot as a PNG image. ::

   png(file="./figure/motifs.png", width = 8, height=ceiling(length(motifs)/3)*3, unit="in", res = 300)
    print(motif_plots)
    dev.off()
  }

Return a list containing the modified Seurat object, the spatial.obj object, and the list of motifs ::

  # return proj_in_tissue
    list(proj_in_tissue = proj_in_tissue, spatial.obj = spatial.obj, motifs = motifs)
  }, error = function(e){
    print(paste0("motif plots skipped. Original error: ", e))
    return(1)
  })
  if(class(out_list) == 'list'){
    proj_in_tissue = out_list$proj_in_tissue
    spatial.obj = out_list$spatial.obj
    motifs = out_list$motifs
    class(proj_in_tissue)
  } else{
    print("out_list not returned")
  }

