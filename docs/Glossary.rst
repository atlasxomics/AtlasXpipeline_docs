Glossary
############

**Data species** refers to the type of organism being analyzed, such as mouse 'mm10', human 'hg38', rat 'rnor6', or a custom reference.

**Number of threads** specifies the number of parallel threads to use when running the analysis, which can improve the efficiency of the analysis by allowing it to run on multiple cores of a computer's processor. An integer must be assigned to num_threads specifying the number of threads to use for ArchR. 

**Tile size** specifies the size of the tiles that will be used to divide the genome into smaller regions for analysis. For tile size we define an integer specifying the bin size of the genome.

**Genome size** specifies the size of the genome being analyzed.

**Minimum TSS and minimum fragments** specify the minimum number of transcription start sites and fragments that must be detected in a region for it to be considered significant. min_Frags refers to the minimum number of mapped ATAC-seq fragments required per tixel when creating arrowFile(s).

**Set seed** is a parameter that can be used to specify a seed value for the random number generator used by the analysis, which can be useful for reproducing results.

**Gene Annotation:** Gene annotation is the process of identifying the locations of genes and other features in a genome, and determining their functions. 

**Arrow Files:** Arrow files are a binary file format used to store columnar data. Arrow files may be used to store the results of gene annotation.

**Maximum Fragments:** Maximum fragments refer to the maximum number of DNA fragments that can be used for gene annotation.

**Add Tile Mat:** The add tile mat function is used to add a tile matrix, which is a matrix of values indicating the presence or absence of a particular feature (such as a gene) in a given region of a genome.

**Add Gene Score Mat:** The add gene score mat function is used to add a gene score matrix, which is a matrix of values indicating the likelihood that a particular region of a genome contains a gene.

**Offset Plus:** Offset plus is a parameter that specifies the number of nucleotides to shift the tile matrix in the positive direction (i.e. to the right) when aligning it with the genome.

**Offset Minus:** Offset minus is a parameter that specifies the number of nucleotides to shift the tile matrix in the negative direction (i.e. to the left) when aligning it with the genome.

**Force:** The force parameter is a logical value that specifies whether to overwrite existing tile or gene score matrices when using the add tile mat or add gene score mat functions.

**Tile Mat Parameters:** Tile mat parameters are a set of parameters that control the behavior of the add tile mat function, such as the offset plus and offset minus values mentioned above. These parameters can be set to customize the alignment of the tile matrix with the genome.



