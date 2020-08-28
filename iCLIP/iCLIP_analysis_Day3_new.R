########################
# load libraries
########################

library(GenomicRanges)
library(rtracklayer)
library(knitr)
library(GenomicFeatures)
library(dplyr)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.NCBI.GRCh38)
library(biomaRt)
library(ggpointdensity)
library(hypeR)
library(eulerr)

## ##############################
## # paths to required files
## ##############################

my_folder = "/Users/kathi/Documents/Teaching/2020_Summer_School_iCLIP/data/"
rbp = "U2AF65"

# file with annotated binding sites
BS_anno_out_file <- paste0(my_folder, "BS_anno_", rbp, ".RData")

# output file in fasta format
fasta_output_file <- paste0("/Users/melinaklostermann/Documents/iCLIP-course/", rbp, "BS_window_40nt.fasta")

# file with annotated binding sites for second RBP
# here: SRSF6
path_RBP2_binding_sites <- paste0(my_folder, "BS_SRSF6.RData")


## ##############################
## # import binding sites
## ##############################

# import annotated binding sites (from Day 2)
load(BS_anno_out_file)

#########################
# prepare binding sites for motif prediction
########################

# use a window of 20 nt on either side
# we take only the top 500 binding sites with highest PureCLIP score

Binding_site_window_for_motif <-
  # make data.frame from GRanges to use arrange function
  as.data.frame(Binding_sites_inside) %>%
  # sort by PureCLIP score
  arrange(desc(score)) %>%
  # subset for top 500 (first 500 rows)
  .[1:500,] %>%
  # turn back into GRanges object
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  # enlarge GRanges by 20 nt on both sites
  + 20

#########################
# get underlying sequences
########################

# load sequences of genome
genome <- BSgenome.Hsapiens.UCSC.hg19 # use GRCH38 for other RBPs

# test getSeq
Binding_sites_inside[1:3]
getSeq(genome, Binding_sites_inside[1:3])

# get sequences for windows around binding sites
Seq_Binding_site_window_for_motif <- getSeq(genome, Binding_site_window_for_motif)
Seq_Binding_site_window_for_motif 

# set name for each sequence in the fasta file
names(Seq_Binding_site_window_for_motif) <- c(seq(1:length(Seq_Binding_site_window_for_motif)))

# export as fasta file
writeXStringSet(Seq_Binding_site_window_for_motif, fasta_output_file)


#########################
# functional annotation of bound genes
########################

# hypeR provides different annotations for several organisms
msigdb_info()

# get pathway info from REACTOME
reactome_geneset <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME")

# load bound genes from Day 2
load(hitgenes_out_file)

# get names of bound genes
bound_gene_names <- hit_genes$gene_name

# test enrichment of genes in Reactome pathways
reactome_hyp <- hypeR(bound_gene_names, reactome_geneset, test="hypergeometric", fdr=0.05)

reactome_hyp

# generate dot plot
reactome_dot <- hyp_dots(reactome_hyp)+
  ggtitle("REACTOME pathways")

reactome_dot


#########################
# overlap binding sites of 2 RBPs
########################

# store binding sites of your RBP in a new variable to avoid mixups
RBP1_binding_sites <- Binding_sites_inside
rm(Binding_sites_inside)

# get file with annotated binding sites for second RBP from your group
# load and rename binding sites of RBP to compare
# here: SRSF6
load(path_RBP2_binding_sites)
RBP2_binding_sites <- Binding_sites_inside

# find overlapping binding sites of both RBPs
overlaps <- subsetByOverlaps(RBP1_binding_sites, RBP2_binding_sites)

# get number of overlapping binding sites
n_overlaps <- NROW(overlaps)

# get total number of binding sites for both RBPs
n_RBP1 <- NROW(RBP1_binding_sites)
n_RBP2 <- NROW(RBP2_binding_sites)

# combine numbers
combined_numbers <-  c(RBP1 = n_RBP1-n_overlaps, RBP2 = n_RBP2-n_overlaps,
                    "RBP1&RBP2" = n_overlaps)

# make plot
fit <- euler(combined_numbers,  shape = "ellipse")
plot(fit, quantities = TRUE)

