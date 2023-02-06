Clustering
###############

**Hierarchical clustering** is a method that uses a tree-like diagram, called a dendrogram, to represent the relationships between different cells. In ArchR, this method is used to group cells based on their similarity in gene expression patterns. The dendrogram shows how the cells are grouped together, with more similar cells being placed closer together in the tree.

**K-means clustering** is a method that involves dividing cells into a specified number of groups, or clusters, based on their similarity. In ArchR, this method is used to group cells based on their gene expression patterns. The user specifies the number of clusters they want to create, and the algorithm will iteratively assign cells to the different clusters in a way that minimizes the distance between cells within a cluster.

**T-SNE (t-distributed stochastic neighbor embedding) clustering** is a method that uses a mathematical technique called dimensionality reduction to project the cells onto a lower-dimensional space, where they can be more easily visualized and clustered. In ArchR, this method is used to group cells based on their gene expression patterns in a way that is particularly useful for visualizing complex data sets with many dimensions.

Clustering Indexes
####################

Calinski-Harabasz index, Dunn index, and Davies-Bouldin index are all measures that can be used to evaluate the quality of clustering results. These measures can be used to compare different clustering methods or different sets of parameters for a particular clustering method, in order to determine which approach produces the most meaningful and interpretable clusters.

**The silhouette width** is a measure of the average distance between a cell and the other cells in its cluster, compared to the average distance between that cell and the cells in the next nearest cluster. A higher silhouette width indicates that the cells in a cluster are more tightly grouped together, and therefore the cluster is more distinct.

**The Calinski-Harabasz index** is a measure of the compactness and separation of clusters. It is calculated as the ratio of the between-cluster variance to the within-cluster variance, where higher values indicate better clustering results.

**The Dunn index** is a measure of the distance between clusters. It is calculated as the minimum distance between any two points in different clusters, divided by the maximum distance between any two points within the same cluster. Higher values of the Dunn index indicate better clustering results, because it means that the clusters are more separated from each other.

**The Davies-Bouldin index** is a measure of the compactness and separation of clusters. It is calculated as the average similarity between each cluster and its most similar neighboring cluster, where lower values indicate better clustering results.

Clustering using Seurat’s FindClusters() function
############################################################

One of the key features of Seurat is its ability to cluster cells based on their gene expression patterns, using a function called FindClusters(). 

The FindClusters() function in Seurat uses a method called graph-based clustering, which involves representing each cell as a point in a high-dimensional space and using a mathematical algorithm to identify clusters of cells that are more similar to each other than they are to cells in other clusters. This method is particularly effective at identifying clusters of cells that have distinct gene expression patterns, even if those patterns are not easily separable using traditional clustering methods.

  projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)
