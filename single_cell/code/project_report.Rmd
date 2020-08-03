---
title: "Project report template"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    code_folding: show
  pdf_document: default
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

<br/>

# Research Analysis Tasks:

## Milestone 1

### Find and download suitable datasets to work
- Consult the Glossary or additional sources for help
- Go to PanglaoDB and look for an appropriate dataset:
- Which species should it be?
- For this tutorial, please use 10X Chromium datasets only
- For this tutorial, please limit the amount of cells to about 10000 cells
- Describe in form of text the rational for this step in your markdown report.


### Load and merge datasets
- Consult the Glossary or additional sources for help
- which file format do we have the data in?
- Why do we need to create a seurat object?
- Describe in form of text the rational for this step in your markdown report.


### Compute QC
- Consult the Glossary or additional sources for help
- Which QC metrics should I calculate?
- What does each QC metric mean?
- Describe in form of text the rational for this step in your markdown report.

### Define appropriate filtering thresholds
- Consult the Glossary or additional sources for help
- Which QC metric needs filtering in your data?
- Which filtering parameters and respective thresholds suit your data?
- Describe in form of text the rational for this step in your markdown report.


<br/>


## Milestone 2

### Normalization and scaling
- Consult the Glossary or additional sources for help
- Why do we need to normalize the data? What exactly are we normalizing?
- Which data covariates could potentially influence the interpreteation of the results?
- Following the question above, are there any covariates that need to be regressed out?
- Are all genes equally important for your analysis? Justify.
- Describe in form of text the rational for this step in your markdown report.


### Data Visualization
- Consult the Glossary or additional sources for help
- Which method would you use to capture most significant information out of the data?
- Following the question above, which parameters would you choose? Why?
- Which method would you choose for visualization of the differences between your cells?
- Following the question above, how the parameters in this method influence your visual representation?
- How some of your QC parameters and datasets influence the separation of your cells?
- Describe in form of text the rational for this step in your markdown report.

<br/>


## Milestone 3

### Dataset integration
- Consult the Glossary or additional sources for help
- Are there batch effects in the data? Is batch correction / dataset integration necessary?
- Following the question above, which parameters allow you to say that batch effects are present?
- Would a simple linear regression be sufficient for removing batch effects on your dataset?
- Which method for batch effect would you choose?
- How could you tell the batch correction procedure worked?
- Describe in form of text the rational for this step in your markdown report.

<br/>


## Milestone 4

### Cell Clustering
- Consult the Glossary or additional sources for help
- What is graph?
- Which kind of graphs is more robust to represent your data?
- Why do we need to cluster our cells?
- Which parameters have you choose when clustering?
- How can you tell which clustering resolution is best?
- Do the clustering reflect the cell separation seen by the visualization method you are using?
- Describe in form of text the rational for this step in your markdown report.

<br/>


## Milestone 5

### Differential expression
- Consult the Glossary or additional sources for help
- Which clustering resolution would you run your differential expression?
- Which test did you choose for differential expression?
- What parameters did you set for computing differential expression? Justify each one.
- Which marker genes can separate each of the cell clusters in your data?
- Which cell types do they represent?
- How would you visualize the lis of differentially expressed genes?
- Describe in form of text the rational for this step in your markdown report.

<br/>


## Milestone 6 (optional)

### Trajectory inference analysis
⁃	Build trajectory embedding (change embedding variables, so you have a continuous path)
⁃	Define lineages
⁃	Visualize lineages
⁃	Perform differential expression on the desired branches
⁃	Report the best gene candidates for each population

