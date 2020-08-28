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
library(UpSetR)


## ##############################
## # paths to required files
## ##############################

my_folder = "/Users/kathi/Documents/Teaching/2020_Summer_School_iCLIP/data/"
rbp = "U2AF65"

## # reimport your Binding sites
path_to_BS_RData <- paste0(my_folder, "Binding_sites_", rbp, ".bed")

# path to single replicates
bw_1_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep1_plus.bw")
bw_1_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep1_minus.bw")
bw_2_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep2_plus.bw")
bw_2_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep2_minus.bw")
bw_3_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep3_plus.bw")
bw_3_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep3_minus.bw")
bw_4_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep4_plus.bw")
bw_4_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep4_minus.bw")
# path to annotation file
path_gtf <- paste0(my_folder, "gencode.v25.annotation.gtf.gz")

# file to save bound genes
hitgenes_out_file <- paste0(my_folder, "bound_genes_", rbp, ".RData")

# file to save annotated binding sites
BS_anno_out_file <- paste0(my_folder, "BS_anno_", rbp, ".RData")


###########################
# import binding sites
###########################

load(path_to_BS_RData)

###########################
# reproducibility of binding sites
###########################

# import crosslink tracks per replicate in Rle format 
sample1.minus.rle <- import.bw( bw_1_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample2.minus.rle <- import.bw( bw_2_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample3.minus.rle <- import.bw( bw_3_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample4.minus.rle <- import.bw( bw_4_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")

sample1.plus.rle <- import.bw( bw_1_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample2.plus.rle <- import.bw( bw_2_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample3.plus.rle <- import.bw( bw_3_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample4.plus.rle <- import.bw( bw_4_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")

# sum up the crosslink events per binding site for each replicate
# plus strand
bs.p = Binding_sites_n_cl[strand(Binding_sites_n_cl) == "+"]
bs.p$clp_rep1 = sample1.plus.rle[bs.p] %>% sum
bs.p$clp_rep2 = sample2.plus.rle[bs.p] %>% sum
bs.p$clp_rep3 = sample3.plus.rle[bs.p] %>% sum
bs.p$clp_rep4 = sample4.plus.rle[bs.p] %>% sum

# minus strand
bs.m = Binding_sites_n_cl[strand(Binding_sites_n_cl) == "-"]
bs.m$clp_rep1 = sample1.minus.rle[bs.m] %>% sum
bs.m$clp_rep2 = sample2.minus.rle[bs.m] %>% sum
bs.m$clp_rep3 = sample3.minus.rle[bs.m] %>% sum
bs.m$clp_rep4 = sample4.minus.rle[bs.m] %>% sum

# combine both strands
bs_cl_samples = c(bs.p, bs.m)

# give every binding an ID
bs_cl_samples$number <- 1:length(bs_cl_samples)

# make list of supported binding sites for the UpSet plot
UpSet_List = list(rep1 = bs_cl_samples[bs_cl_samples$clp_rep1>0]$number,
                rep2 = bs_cl_samples[bs_cl_samples$clp_rep2>0]$number,
                rep3 = bs_cl_samples[bs_cl_samples$clp_rep3>0]$number,
                rep4 = bs_cl_samples[bs_cl_samples$clp_rep4>0]$number)

UpSetR::upset(UpSetR::fromList(UpSet_List), order.by = c("freq", "degree"), nsets = 4)


###########################
# scatter plots between replicates
###########################

# scatter plot for replicates 1 and 2
head(bs_cl_samples)
bs_cl_samples_df <- as.data.frame(bs_cl_samples)

# make plot
ggplot(bs_cl_samples_df, aes(x=clp_rep1, y=clp_rep2))+
  geom_point()

# log scale  
ggplot(bs_cl_samples_df, aes(x=log2(clp_rep1), y=log2(clp_rep2)))+
  geom_point()

# point density (this is nice, but very compute-intensive!)
#ggplot(bs_cl_samples_df, aes(x=log2(clp_rep1), y=log2(clp_rep2)))+
#  geom_pointdensity()

# calculate Pearson correlation
ggplot(bs_cl_samples_df, aes(x=log2(clp_rep1), y=log2(clp_rep2)))+
  geom_point()+
  stat_cor()

###########################
# filter for reproducible binding sites
###########################

bs_cl_samples_df <- bs_cl_samples_df %>% rowwise %>%
  mutate(n_repro_reps = sum(clp_rep1>0, 
                             clp_rep2>0,
                             clp_rep3>0,
                             clp_rep4>0))

head(bs_cl_samples_df)

# how many binding sites are supported by less than 3 replicates?
nrow(bs_cl_samples_df[bs_cl_samples_df$n_repro_reps<3,])

# filter for 3 or more replicates
bs_cl_samples <- bs_cl_samples[bs_cl_samples_df$n_repro_reps>=3]

# use these as binding sites in the following analyses
Binding_sites <- bs_cl_samples


###########################
# import gene annotation
###########################

# load annotation
gtf <- import(path_gtf)
gtf

# my favorite gene is PURA, the Ensembl gene_id is ENSG00000185129
# I can subset the GRanges objects by choosing all rows that have PURA in the gene_name column
anno_favorite_gene <- gtf[gtf$gene_name=="PURA"]
anno_favorite_gene

# with the table() function I can look at how many annotations of which type are there for PURA
table(anno_favorite_gene$type)

# GRanges offers a function that counts overlaps 
countOverlaps(anno_favorite_gene)

# this tells us how many entries overlap with each row
# with the first row of anno_favorite_gene (this is the row of the gene) overlap 14 other entries (all)
# with the second row (one transcript) overlap 12 entries and so on

###########################
# filter gene annotation
###########################

# select only genes from annotations
anno_genes <- gtf[gtf$type=="gene"]

# are ther overlaps?
head(countOverlaps(anno_genes), 20)

# yes many have overlaps eg. the gene in the 2. and the 3. row of the anno_genes granges
anno_genes[2:3]

#  first we filter the annotation for standard chromosomes
gtf <-keepStandardChromosomes(gtf, pruning.mode = "coarse")
  
# select genes from the annotation
gtf_genes <- gtf[gtf$type=="gene"]
# keep all genes with level 1 or 2
gtf_genes_GL12<- gtf_genes[gtf_genes$level <= 2]
  
# keep gene with level 3 only if they do not overlap
gtf_GL3 <- subsetByOverlaps(gtf_genes[gtf_genes$level==3], gtf_genes_GL12, type = "any", invert = T)
genes <- c(gtf_genes_GL12, gtf_GL3)  


######################
# plot bound genes
#####################

# subset genes with binding sites
hit_genes <- subsetByOverlaps(genes, Binding_sites)
hit_genes

# save bounds genes for GO analysis on Day 3
save(hit_genes, file = hitgenes_out_file)


# how many genes are hit?
NROW(hit_genes)

# make a data frame for ggplot that contains the gene_types
hit_genetypes_df <- data.frame(gene_types=hit_genes$gene_type)

# plot gene_types
ggplot(hit_genetypes_df , aes(x = gene_types)) + 
  geom_bar(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # turns axis label 90 degree
  ggtitle("Bound genes") + # give plot a title
  xlab("labels") + # change label of x-axis
  ylab("counts [log10]") # and of y-axis

# to see all bars, we can use a log scale on the y-axis
ggplot(hit_genetypes_df , aes(x = gene_types)) + 
  geom_bar(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Bound genes (log)") + 
  scale_y_log10() + 
  xlab("labels") + 
  ylab("counts [log10]")


###########################
# extract transcript regions
###########################

# make txDb (GenomicFeatures package)
anno_txDb <- makeTxDbFromGRanges(gtf)
anno_txDb

# extract different transcript regions
# get introns
introns <- intronsByTranscript(anno_txDb)
introns #this returns GRangesList, but we want to have GRanges again as output
introns <- unlist(introns) # this converts GRangesList into GRanges again
introns

# 3'UTR
utrs3 <- unlist(threeUTRsByTranscript(anno_txDb))

# 5'UTR
utrs5 <- unlist(fiveUTRsByTranscript(anno_txDb))

# CDS
cds <- cds(anno_txDb)


######################### 
### conflicting annotations - problem visualization
######################## 

# count overlap of each binding site within each part of the gene
cds_BS <- countOverlaps(Binding_sites, cds)
introns_BS <- countOverlaps(Binding_sites, introns)
utrs3_BS <- countOverlaps(Binding_sites, utrs3)
utrs5_BS <- countOverlaps(Binding_sites, utrs5)

count.df <- data.frame(cds = cds_BS, intron = introns_BS, utr3 = utrs3_BS, utr5 = utrs5_BS)
head(count.df)

### plot number of different transcript annotations at same position
n_overlaps = apply(count.df, 1, function(x) length(x[x != 0])) # count how many regions overlap with each binding site
head(df)

# data.frame for ggplot
df <- data.frame(n_overlaps = n_overlaps)

# plot
ggplot(df, aes(x = n_overlaps)) + 
  geom_bar(stat = "count") + 
  ggtitle("Different transcript region overlaps") + 
  xlab("number of different annotations") +
  ylab("counts")

######################### 
### resolve conflicting annotations
######################## 

# Set hierarchical rule for ties
rule = c("utr3", "intron", "cds", "utr5") # change this according to your rule

# apply majority vote
count.df.reg = count.df[, rule] %>% # this orders the columns in the order of our rule
 mutate(., outside = case_when(rowSums(.) == 0 ~ TRUE,
                               TRUE~FALSE))%>% 
rowwise() %>%
  mutate(region = which.max(c(utr3,intron, cds, utr5))) %>%
  ungroup() %>%
  mutate( region = case_when( region == 1 ~ rule[1],
                              region == 2 ~ rule[2],
                              region == 3 ~ rule[3],
                              region == 4 ~ rule[4])
          )
  
head(count.df.reg)

Binding_sites$outside = count.df.reg$outside
Binding_sites$region = count.df.reg$region

Binding_sites_inside = Binding_sites[Binding_sites$outside==FALSE]

# make plot
df = data.frame(region = Binding_sites_inside$region)
ggplot(df, aes(x = region)) +
  geom_bar() +
  xlab("region") +
  ylab("count")


######################### 
### save annotated binding sites
######################## 

save(Binding_sites_inside, file = BS_anno_out_file)

