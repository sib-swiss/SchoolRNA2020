---
title: "Glossary of terms"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    code_folding: show
  pdf_document: default
editor_options: 
  chunk_output_type: console
---

<style>

body,p,h1,h2,h3,h4,h5,h6 {
  text-align: justify;
  color: #424949;
  font-weight: 400;
  font-size: 14px;
  line-height: 1.55;}

h1 {  font-size: 1.4em;font-weight: bold; }
h2 {  font-size: 1.2em;font-weight: bold; }
h3 {  font-size: 1.1em;font-weight: bold; }
h4 {  font-size: 1.0em;font-weight: bold; }

pre {
  font-family: "Courier";
  background-color: #F7F7F7;
  border: 1px solid #CCCCCC;
  font-size: 14px;}

code {
  font-family: "Courier New";
  font-size: 14px;
  word-break: break-all;
  font-weight: bold;}

</style>



# Quality control
***

<br/>

## Total no. of features

A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. Please note that those values are highly dependent on the library preparation method used.

<br/>

## Number of features per cell

<br/>

## % Mitochondrial genes

Having the data in a suitable format, we can start calculating some quality metrics. We can for example calculate the percentage of mitocondrial and ribosomal genes per cell and add to the metadata. This will be helpfull to visualize them across different metadata parameteres (i.e. datasetID and chemistry version). There are several ways of doing this, and here manually calculate the proportion of mitochondrial reads and add to the metadata table.

Citing from “Simple Single Cell” workflows (Lun, McCarthy & Marioni, 2017): “High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.”


```r
# Calculating % mitochondrial genes
SeuratObject <- PercentageFeatureSet(SeuratObject, 
                                     pattern = "^MT-",
                                     col.name = "percent_mito")

# Filter
selected_mito <- WhichCells(SeuratObject,
                            expression = percent_mito < 0.25)
SeuratObject <- subset(SeuratObject,
                       cells = selected_mito)
```

<br/>

## % Ribossomal genes


```r
# Calculating % ribossomal genes
SeuratObject <- PercentageFeatureSet(SeuratObject, 
                                     pattern = "^RP[SL]", 
                                     col.name = "percent_ribo")

# Filter
selected_ribo <- WhichCells(SeuratObject, 
                            expression = percent_ribo > 0.05)
SeuratObject <- subset(SeuratObject, 
                       cells = selected_ribo)
```

<br/>

## Removal of genes

In RNA-sequencing, genes can be categorized into different groups depending on their RNA biotype. For example, "coding", "non-coding", "VDJ region genes" are "small interefering RNA" common gene biotypes. Besides, having information about chromossomal location might be usefull to identify bacth effects driven by sex chromossomes.

Depending on the desired type of analysis, some gene categories can be filtered out if not of interest. For single cell specifically, cell libraries are usually constructed using poly-A enrichment and therefore enriching for "protein-coding proteins", which usually contitutes around 80-90% of all available genes.

How to run it:


```r
mart = useMart("ensembl", 
               dataset = paste0("mmusculus_gene_ensembl"),
               host="jul2019.archive.ensembl.org")
annot <- getBM(c("external_gene_name",
                 "gene_biotype",
                 "chromosome_name"),
               mart = mart)

sel <- annot[match(rownames(DATA@assays[["RNA"]]@counts) , annot[,1]),2] == "protein_coding"
genes_use <- rownames(DATA@assays[["RNA"]]@counts)[sel]
genes_use <- as.character(na.omit(genes_use))
DATA@assays[["RNA"]]@counts <- DATA@assays[["RNA"]]@counts[genes_use,]
```

## Cell cycle

We here perform cell cycle scoring. To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list. Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.

How to run it:


```r
SeuratObject <- CellCycleScoring(object = SeuratObject,
                                 g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
```

Plot:


```r
VlnPlot(SeuratObject,
        features = c("S.Score","G2M.Score"),
        group.by= "orig.ident",
        ncol = 4,
        pt.size = .1)
```


<br/>

<br/>

# Normalization and Regression
***

## Normalization



How to run it:


```r
SeuratObject <- NormalizeData(object = SeuratObject,
                              scale.factor = 10000,
                              normalization.method = "LogNormalize")
```


## Scaling and Centering (linear)

Since each gene has a different expression level, it means that genes with higher expression values will naturally have higher variation that will be captured by downstream methods. This means that we need to somehow give each gene a similar weight beforehand (see below). A common practice is to center and scale each gene before performing PCA. This exact scaling is called Z-score normalization it is very useful for PCA, clustering and plotting heatmaps.

Additionally, we can use regression to remove any unwanted sources of variation from the dataset, such as cell cycle, sequencing depth, percent mitocondria. This is achieved by doing a generalized linear regression (GLM) using these parameters as covariates in the model. Then the residuals of the model are taken as the “regressed data”. Although perhaps not in the best way, batch effect regression can also be done here.

How to run it:


```r
SeuratObject <- ScaleData(object = SeuratObject,
                          vars.to.regress = c("nUMI","mito.percent","nFeatures"),
                          model.use = "linear",
                          do.scale = T,
                          do.center = T)
```


## Scaling and Centering (poisson)

Since the procedure above assumes a log-linear data distribution, it may be the case that it does not regress the variation correctly, as RNA-seq data (including single cell) relates more closely to a negative bionomial distribution. An alternative variation of the procedure above can also be run on the raw UMI count data but using a "poisson" or "negative binomial" distribution instead. This is performing a gene-wise GLM regression using a poisson model.

How to run it:


```r
SeuratObject <- ScaleData(object = SeuratObject,
                          vars.to.regress = c("nUMI","mito.percent","nFeatures"),
                          model.use = "poisson",
                          do.scale = T,
                          do.center = T)
```

## SCtransform

Scaling and centering assuming a poisson distribution might in some cases overfit the data, see above. One can overcome this by pooling information across genes with similar abundances in order to obtain more stable parameter estimates to be used as gene weights in the regression model. This is called "scTransform" and, in simple terms, is performing a gene-wise GLM regression using a contrained negative binomial model.

How to run it:


```r
SeuratObject <- SCTransform( object = SeuratObject,
                             assay=opt$assay,
                             new.assay.name = "sctransform",
                             do.center=T )
```

# Feature selection
***

An important step in many big-data analysis tasks is to identify features (genes, transcripts, proteins, metabolites, etc) that are actually very variable between the samples being looked at.

For example. Imagine that you have a dataset known to contain different types of epithelial cells, and you use either 1) only genes that are expressed and shared across all epithelial cells at about the same level, 2) only genes that are not detected in epithelial cells, 3) only genes which expression differ greatly across epithelial cells or 4) using all genes. Which of these 4 gene lists can best distinguish the epithelial subtypes in this dataset?

As you could now imagine, using only genes which expression differ greatly across epithelial cells is the best case scenario, followed by using al genes. Therefore, using only genes that are expressed and shared across all epithelial cells at about the same level or only genes that are not detected in epithelial cells do not contain sufficient information to distinguish the epithelial subtypes.

However, since in single-cell we usually do not know the epithelial suptypes the cells before hand (since this is what we want to discover), we need another method to acomplish this task. In general terms, a common approach is to order genes by their overal variance across samples. This is because genes with higher variance will also likely be the ones that can separate the cells the best.

Since genes with higher expression level usually also have naturally higher variation, the gene variation is then normalized by the log  mean expression of each gene (see plot). 

How to run it:


```r
SeuratObject <- FindVariableFeatures(
  object = SeuratObject,
  nfeatures = 2000,
  selection.method = "vst",
  verbose = FALSE,
  assay = "RNA",
  dispersion.function = FastLogVMR,
  mean.function = FastExpMean)
```

Variable gene plot:


```r
top20 <- head(VariableFeatures(alldata), 20)
LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
```

# Intro to Graphs
***

## KNN

## SNN


How to run it:

```r
SeuratObject <- FindNeighbors(SeuratObject,
                              assay = "RNA",
                              reduction = "pca" 
                              dims = 1:50,
                              graph.name="SNN",
                              prune.SNN = 1/15,
                              k.param = 20,
                              force.recalc = T)
```

<br/>

<br/>

# Dimensionality reduction
***

<br/>

## PCA

Principal Component Analysis (PCA) is defined as an orthogonal **linear** transformation that transforms the data to a new coordinate system such that **the greatest variance by some scalar projection of the data comes to lie on the first coordinate** (called the first principal component), the second greatest variance on the second coordinate, and so on. […] Often, its operation can be thought of as revealing the internal structure of the data in a way that best explains the variance in the data. […] This is done by **using only the first few principal components** so that the dimensionality of the transformed data is reduced.

<div style="text-align: right"> Adapted from [Wikipedia](https://en.wikipedia.org/wiki/Principal_component_analysis) </div>

How to run it:

```r
SeuratObject <- RunPCA(object = SeuratObject,
                       assay = "RNA",
                       npcs = 100,
                       verbose = FALSE )
```

<br/> 

<br/> 

## tSNE

<div style="text-align: right"> [Maaten, Hilton (2008) J of Machine Learning Research](http://jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf) </div>


t-distributed stochastic neighborhood embedding (tSNE) is a **nonlinear** dimensionality reduction technique well-suited for embedding high-dimensional data for **visualization** in a low-dimensional space of two or three dimensions. Specifically, it models each high-dimensional object by a two- or three-dimensional point in such a way that **similar objects are modeled by nearby points** and dissimilar objects are modeled by distant points with high probability. […] t-SNE has been used for visualization in a wide range of applications, including […] bioinformatics […]. While t-SNE plots often seem to display clusters, the **visual clusters can be influenced strongly by the chosen parameterization** and therefore a good understanding of the parameters for t-SNE is necessary. 

<div style="text-align: right"> Adapted from [Wikipedia](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) </div>


Usefull links:

* [How to Use t-SNE Effectively](https://distill.pub/2016/misread-tsne/)

How to run it:

```r
SeuratObject <- RunTSNE(object = SeuratObject,
                        reduction = "pca",
                        perplexity=30,
                        max_iter=1000,
                        theta=0.5,
                        eta=200,
                        exaggeration_factor=12,
                        dims.use = 1:50,
                        verbose = T,
                        num_threads=0)
```

<br/> 

<br/> 

## UMAP

Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general **non-linear** dimension reduction […]. 

<div style="text-align: right"> [umap-learn documentation](https://umap-learn.readthedocs.io/en/latest/) </div>

The result is a practical scalable algorithm that applies to real world data. The UMAP algorithm is competitive with t-SNE for **visualization** quality, and arguably preserves **more of the global structure** with superior run time performance. Furthermore, UMAP has no computational restrictions on embedding dimension, making it viable as a general purpose dimension reduction technique for machine learning. 

<div style="text-align: right"> [UMAP Arxiv paper](https://arxiv.org/pdf/1802.03426.pdf) </div>

How to run it:

```r
SeuratObject <- RunUMAP(object = SeuratObject,
                        reduction = "pca",
                        dims = 1:top_PCs,
                        n.components = 2,
                        n.neighbors = 20,
                        spread = .3,
                        repulsion.strength = 1,
                        min.dist= .001,
                        verbose = T,
                        num_threads=0,
                        n.epochs = 200,
                        metric = "euclidean",
                        seed.use = 42,
                        reduction.name="umap")
```

<br/> 

<br/> 

## DM

Diffusion maps (DM) is a dimensionality reduction [...] which computes a family of embeddings of a data set into Euclidean space (often low-dimensional) whose coordinates can be computed from the eigenvectors and eigenvalues of a diffusion operator on the data. The Euclidean distance between points in **the embedded space is equal to the "diffusion distance" between probability distributions** centered at those points. Different from linear dimensionality reduction methods such as principal component analysis (PCA) and multi-dimensional scaling (MDS), diffusion maps is part of the family of **nonlinear** dimensionality reduction methods which focus on discovering the underlying manifold that the data has been sampled from. [...] The basic observation is that **if we take a random walk on the data, walking to a nearby data-point is more likely than walking to another that is far away**.

<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/Diffusion_map) </div>


[Diffusion Maps paper](https://www.pnas.org/content/pnas/102/21/7426.full.pdf)

How to run it:

```r
# Load additional libraries
library(destiny)

#Run diffusion maps using the destiny package 
dm <- DiffusionMap( data = SeuratObject@reductions[["pca"]]@cell.embeddings[ , 1:50],
                    k = 20,
                    n_eigs = 20)

#Fix the cell names in the DM embedding
rownames(dm@eigenvectors) <- colnames(SeuratObject)

#Add the DM embbedding to the SeuratObject
SeuratObject@reductions[["dm"]] <- CreateDimReducObject(embeddings = dm@eigenvectors,
                                                        key = "DC_",
                                                        assay = "RNA")
```

<br/> 

<br/> 


## ICA

Independent Component Analysis (ICA) is a computational method for separating a multivariate signal into additive subcomponents. This is done by assuming that the subcomponents are non-Gaussian signals and that they are statistically independent from each other. ICA is a special case of blind source separation.

<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/Independent_component_analysis) </div>

How to run it:

```r
SeuratObject <- RunICA(object = SeuratObject,
                       assay = "pca",
                       nics = 20,
                       reduction.name = "ica")
```

<br/> 

<br/> 


# Dataset integration
***

<br/> 

## MNN

<div style="text-align: right"> [Haghverdi et al (2018) *Nat Biotechnology*](https://www.nature.com/articles/nbt.4091) </div>


```r
# Load additional libraries
library(SeuratWrappers)

SeuratObject.list <- SplitObject(SeuratObject, split.by = "BATCH")
SeuratObject <- RunFastMNN(object.list = SeuratObject.list,
                           assay = "RNA",
                           features = 2000,
                           reduction.name = "mnn")

# Free memory from working environment
rm( c( SeuratObject.list ) )
gc(verbose = FALSE)
```

## CCA

<div style="text-align: right"> [Stuart et al (2019) *Cell*](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) </div>


```r
SeuratObject.list <- SplitObject(SeuratObject, split.by = "BATCH")

for (i in 1:length(SeuratObject.list)) {
    SeuratObject.list[[i]] <- NormalizeData(SeuratObject.list[[i]], verbose = FALSE)
    SeuratObject.list[[i]] <- FindVariableFeatures(SeuratObject.list[[i]],
                                                   selection.method = "vst",
                                                   nfeatures = 2000,
                                                   verbose = FALSE)
}

hvgs_per_dataset <- lapply(SeuratObject.list, function(x) { x@assays$RNA@var.features })

venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn = 1,cexil = 1,lwd=1,col="white",frame=F,borders = NA)

SeuratObject.anchors <- FindIntegrationAnchors(object.list = SeuratObject.list, dims = 1:30)

SeuratObject.int <- IntegrateData(anchorset = SeuratObject.anchors,
                                  dims = 1:30,
                                  new.assay.name = "cca")
```

<br/> 

<br/> 

## LIGER

<div style="text-align: right"> [Welch et al (2019) *Cell*](https://www.cell.com/cell/pdf/S0092-8674(19)30504-5.pdf) </div>


```r
# Load additional libraries
library(conos)
library(SeuratWrappers)

# Split the data per batch to be corrected
SeuratObject.list <- SplitObject(SeuratObject, split.by = "Method")
for (i in 1:length(SeuratObject.list)) {
    SeuratObject.list[[i]] <- 
      NormalizeData(SeuratObject.list[[i]]) %>%
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA(verbose = FALSE)
}

# Create a Conos object
SeuratObject.con <- Conos$new(SeuratObject.list)

# Build a joint graph across datasets and find shared communities
SeuratObject.con$buildGraph(k = 15, 
                            k.self = 5,
                            space = "PCA",
                            ncomps = 30,
                            n.odgenes = 2000,
                            matching.method = "mNN",
                            metric = "angular",
                            score.component.variance = TRUE,
                            verbose = TRUE)
SeuratObject.con$findCommunities()

# Create a Joint embedding and conver it back to Seurat Object
SeuratObject.con$embedGraph()
SeuratObject <- as.Seurat(SeuratObject.con)

# Free memory from working environment
rm( c( SeuratObject.con, SeuratObject.list ) )
gc(verbose = FALSE)
```

<br/> 

<br/> 

## Conos

<div style="text-align: right"> [Barkas et al (2019) *Nat Methods*](https://www.nature.com/articles/s41592-019-0466-z) </div>


```r
# Load additional libraries
library(conos)
library(SeuratWrappers)

# Split the data per batch to be corrected
SeuratObject.list <- SplitObject(SeuratObject, split.by = "Method")
for (i in 1:length(SeuratObject.list)) {
    SeuratObject.list[[i]] <- 
      NormalizeData(SeuratObject.list[[i]]) %>%
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA(verbose = FALSE)
}

# Create a Conos object
SeuratObject.con <- Conos$new(SeuratObject.list)

# Build a joint graph across datasets and find shared communities
SeuratObject.con$buildGraph(k = 15, 
                            k.self = 5,
                            space = "PCA",
                            ncomps = 30,
                            n.odgenes = 2000,
                            matching.method = "mNN",
                            metric = "angular",
                            score.component.variance = TRUE,
                            verbose = TRUE)
SeuratObject.con$findCommunities()

# Create a Joint embedding and conver it back to Seurat Object
SeuratObject.con$embedGraph()
SeuratObject <- as.Seurat(SeuratObject.con)

# Free memory from working environment
rm( c( SeuratObject.con, SeuratObject.list ) )
gc(verbose = FALSE)
```

<br/> 

<br/> 


## Harmony

<div style="text-align: right"> [Korsunsky et al (2019) *Nat Mathods*](https://www.nature.com/articles/s41592-019-0619-0) </div>



```r
# Load additional libraries
library(harmony)
library(SeuratWrappers)

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = "Method")
```

<br/> 

<br/> 


# Clustering
***

<br/> 

## Louvain

The Louvain method for community detection is a method to extract communities from large networks created by Blondel et al. from the University of Louvain. The method is a greedy optimization method that appears to run in time $O(n.log^2n)$ in the number of nodes in the network.The value to be optimized is modularity, defined as a value in the range that measures the density of links inside communities compared to links between communities. Optimizing this value theoretically results in **the best possible grouping of the nodes of a given network**, however going through all possible iterations of the nodes into groups is impractical so heuristic algorithms are used.

<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/Independent_component_analysis) </div>


[Louvain Paper](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/pdf)

How to run it:

```r
SeuratObject <- RunICA(object = SeuratObject,
                       assay = "pca",
                       nics = 20,
                       reduction.name = "ica")
```

<br/> 

<br/> 

## Leiden

Leiden algorithm is applied iteratively, it converges to a partition in which all subsets of all communities are locally optimally assigned. Furthermore, by relying on a fast local move approach, the Leiden algorithm runs faster than the Louvain algorithm. The Leiden algorithm consists of three phases: (1) local moving of nodes, (2) refinement of the partition and (3) aggregation of the network based on the refined partition, using the non-refined partition to create an initial partition for the aggregate network.

<div style="text-align: right"> [Leiden Paper](https://www.nature.com/articles/s41598-019-41695-z.pdf) </div>

<br/> 

<br/> 


## Hierachical clustering

Hierachical clustering (HC) is a method of cluster analysis which **seeks to build a hierarchy of clusters**. Strategies for hierarchical clustering generally fall into two types: Agglomerative or Divisive. [...] In general, the merges and splits are determined in a greedy manner. The results of hierarchical clustering are **usually presented in a dendrogram**. [...] The standard algorithm for hierarchical agglomerative clustering (HAC) has a time complexity of $O(n^3)$ and requires $O(n^2)$ memory, which makes it **too slow for even medium data sets**. In order to decide which clusters should be combined (for agglomerative), or where a cluster should be split (for divisive), **a measure of dissimilarity between sets of observations** is required. In most methods of hierarchical clustering, this is achieved by use of an appropriate metric (a measure of distance between pairs of observations), and **a linkage criterion** which specifies the dissimilarity of sets as a function of the pairwise distances of observations in the sets.

<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/Hierarchical_clustering) </div>

[HC for networks](https://en.wikipedia.org/wiki/Hierarchical_clustering_of_networks)

<br/> 

<br/> 

## K-means

k-means clustering is a method of vector quantization, originally from signal processing, that aims to partition $n$ observations into $k$ clusters in which **each observation belongs to the cluster with the nearest mean** (cluster centers or cluster centroid), serving as a prototype of the cluster. [...] The algorithm does not guarantee convergence to the global optimum. The result may depend on the initial clusters. As the algorithm is usually fast, **it is common to run it multiple times** with different starting conditions. [...] k-means clustering tends to find clusters of **comparable spatial extent (all with same size)**, while the expectation-maximization mechanism allows clusters to have different shapes.
<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/K-means_clustering) </div>

<br/> 

<br/> 

## Affinity propagation

In statistics and data mining, affinity propagation (AP) is a clustering algorithm based on the concept of "message passing" between data points. Unlike clustering algorithms such as k-means or k-medoids, affinity propagation **does not require the number of clusters to be determined** or estimated before running the algorithm. Similar to k-medoids, affinity propagation finds "exemplars," members of the input set that are representative of clusters. [...] Iterations are performed until either the cluster boundaries remain unchanged over a number of iterations, or some predetermined number (of interations) is reached.

<div style="text-align: right"> [Wikipedia](https://en.wikipedia.org/wiki/Affinity_propagation) </div>

<br/> 

<br/> 


# Trajectory inference

## Affinity propagation


## Affinity propagation




