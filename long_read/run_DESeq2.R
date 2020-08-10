cts <- read.delim("~/flair_output/counts_matrix.tsv", row.names = 1)
coldata <- read.delim("~/reads_manifest.tsv", header = F)
coldata$V1 <- gsub("dorsolateral_prefrontal", "dlp", coldata$V1)
tiss_subj <- do.call(rbind, strsplit(coldata$V1, '-'))
coldata <- data.frame(tissue = factor(tiss_subj[,1]), subject = factor(tiss_subj[,2]), row.names = coldata$V1)

cts["ENST00000399641_ENSG00000151067",]

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ tissue + subject)
dds <- DESeq(dds)

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds)
rld <- rlog(dds, blind=FALSE)

library("vsn")
meanSdPlot(assay(ntd))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("subject","tissue")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

df <- as.data.frame(colData(dds)[,c("tissue","subject")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(ntd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tissue, vsd$subject, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

DESeq2::plotPCA(vsd, intgroup=c("tissue"))
DESeq2::plotPCA(vsd, intgroup=c("subject"))

coldata2 <- coldata
coldata2$tissue <- as.character(coldata2$tissue)
ewc <- endsWith(coldata2$tissue, "cortex")
coldata2$tissue[ewc] <- "cortex"

dds2 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata2,
                              design = ~ tissue)
dds2 <- DESeq(dds2)

res <- results(dds2, contrast = c("tissue", "cerebellum", "cortex"))

plotMA(res, ylim=c(-3,3))

plotCounts(dds2, gene=which.min(res$padj), intgroup="tissue")

maxg <- res@rownames[which.min(res$padj)]
cts[maxg,]
rel_counts <- t(apply(cts, 1, function(x){
  x/colSums(cts)
}))

relc_max <- rel_counts[maxg,]
relc_high <- rel_counts["ENST00000399641_ENSG00000151067",]

boxplot(relc_max ~ coldata2$tissue)
boxplot(relc_high ~ coldata2$tissue)
