# Exercises day 3

In the below exercises we will analyse the differential isoform expression. To do this, we will use the R package DESeq2. For more information, read the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

In this markdown file, we'll showcase a basic analysis as can be done by DESeq2. However, if there's time left, try to go a bit deeper and try for example different ways of visualisation, normalisation or experimental design.

## 3.1 Run DESeq2
><img border="0" src="https://www.svgrepo.com/show/14756/person-silhouette.svg" width="30" height="30"> 1 hour

Load the data into R and modify it to fit DESeq2:

```r
# read count data and metadata (coldata)
cts <- read.delim("~/flair_output/counts_matrix.tsv", row.names = 1)
coldata <- read.delim("/data/reads/lrrnaseq/batch_combined/reads_manifest.tsv", header = F)

# make sample names a bit shorter:
coldata$V1 <- gsub("dorsolateral_prefrontal", "dlp", coldata$V1)

# modify coldata to fit DESeq2
tiss_subj <- do.call(rbind, strsplit(coldata$V1, '-'))
coldata <- data.frame(tissue = factor(tiss_subj[,1]), subject = factor(tiss_subj[,2]), row.names = coldata$V1)

head(cts[,1:4])
coldata

```

Make a DESeq2 object called `dds`:

```r
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ tissue + subject)
```

Normalize the counts for visualisation:

```r
ntd <- DESeq2::normTransform(dds)
```

Let's start with a PCA:

```r
DESeq2::plotPCA(ntd, intgroup=c("tissue"))
```

It seems that cortex samples cluster, and cerebellum is most different from the rest. Let's change the experimental design:

```r
coldata2 <- coldata
coldata2$tissue <- as.character(coldata2$tissue)
ewc <- endsWith(coldata2$tissue, "cortex")
coldata2$tissue[ewc] <- "cortex"
coldata2$tissue <- factor(coldata2$tissue)
coldata2
```

And make a new DESeq2 object

```r
dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata2,
                              design = ~ tissue)
```

Let's compare all cortex samples with cerebellum samples:

```r
dds2 <- DESeq2::DESeq(dds2)

res <- DESeq2::results(dds2, contrast = c("tissue", "cerebellum", "cortex"))

res
```

And let's see which isoform is most significantly different between the two:

```r
most_significant_isoform <- res@rownames[which.min(res$padj)]
most_significant_isoform

```

And visualise the difference between the tissues in expression:

```r
DESeq2::plotCounts(dds2, gene=most_significant_isoform, intgroup="tissue")
```

## 3.2 Biological meaning
><img border="0" src="https://www.svgrepo.com/show/220819/group-team.svg" width="30" height="30"> 2 hours

Try to answer one or more of the following questions:
* Which isoform is most differentially expressed between cerebellum and cortex? Is it a new isoform? Did you find the same as in the publication?
* Is there a difference between the relative number of productive expressed isoforms between the tissues?
* How are the samples correlated? Is there a hierarchical clustering as expected?

Summarize your findings of the entire project in a 20 minute presentation to show to the other participants.

### Go to:
* [Main page](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq)
