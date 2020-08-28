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
library(ggpointdensity)

## ##############################
## # paths to required files
## ##############################

my_folder = "/Users/kathi/Documents/Teaching/2020_Summer_School_iCLIP/data/"
rbp = "U2AF65"

# path to single replicates
bw_1_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep1_plus.bw")
bw_1_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep1_minus.bw")
bw_2_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep2_plus.bw")
bw_2_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep2_minus.bw")
bw_3_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep3_plus.bw")
bw_3_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep3_minus.bw")
bw_4_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep4_plus.bw")
bw_4_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_rep4_minus.bw")

# path to merged bigWig files
bw_plus_path <- paste0(my_folder, "crosslinks_",rbp,"_all_plus.bw")
bw_minus_path <- paste0(my_folder, "crosslinks_",rbp,"_all_minus.bw")

# path to PureCLIP calls
pureclip_path <- paste0(my_folder, "peakcalling_",rbp,"_sites.bed")

# output file for intermediate binding sites
bs_intermediate_out_file <- paste0(my_folder, "binding_sites_intermediate_",rbp,".bed")

# output file for resized binding sites
bs_out_file <- paste0(my_folder, "Binding_sites_", rbp, ".bed")

########################
# import bigWig files (crosslink events)
########################

# import replicate 1
bw_1_plus <- import.bw(bw_1_plus_path)
bw_1_minus <- import.bw(bw_1_minus_path)

# show imported files
bw_1_plus
bw_1_minus

# replicate 2
bw_2_plus <- import.bw(bw_2_plus_path)
bw_2_minus <- import.bw(bw_2_minus_path)

bw_2_plus
bw_2_minus

# replicate 3
bw_3_plus <- import.bw(bw_3_plus_path)
bw_3_minus <- import.bw(bw_3_minus_path)

bw_3_plus
bw_3_minus

# replicate 4
bw_4_plus <- import.bw(bw_4_plus_path)
bw_4_minus <- import.bw(bw_4_minus_path)

bw_4_plus
bw_4_minus

# add strand information to the GRanges object
strand(bw_1_plus) <- "+"
strand(bw_1_minus) <- "-"
# and have a look again:
bw_1_plus
bw_1_minus
# now the strands are annotated correctly with "+" or "-"

# add_strand_others
strand(bw_2_plus) <- "+"
strand(bw_2_minus) <- "-"

strand(bw_3_plus) <- "+"
strand(bw_3_minus) <- "-"

strand(bw_4_plus) <- "+"
strand(bw_4_minus) <- "-"


#######################
# remove additional scaffolds
#######################

# look at seqnames 
seqnames(bw_1_plus)

# replicate 1
bw_1_plus <- keepStandardChromosomes(bw_1_plus, pruning.mode = "coarse")
bw_1_minus <- keepStandardChromosomes(bw_1_minus, pruning.mode = "coarse")

bw_1_plus
bw_1_minus

# other replicates
bw_2_plus <- keepStandardChromosomes(bw_2_plus, pruning.mode = "coarse")
bw_2_minus <- keepStandardChromosomes(bw_2_minus, pruning.mode = "coarse")

bw_3_plus <- keepStandardChromosomes(bw_3_plus, pruning.mode = "coarse")
bw_3_minus <- keepStandardChromosomes(bw_3_minus, pruning.mode = "coarse")

bw_4_plus <- keepStandardChromosomes(bw_4_plus, pruning.mode = "coarse")
bw_4_minus <- keepStandardChromosomes(bw_4_minus, pruning.mode = "coarse")

#######################
# collect some stats about data set
#######################

# number of crosslink sites
NROW(bw_1_plus)

# number of crosslink events
sum(bw_1_plus$score)

# average crosslink events per site
sum(bw_1_plus$score)/NROW(bw_1_plus)

# make a table of crosslink sites and events of all replicates
crosslinks_overview_table <- data.frame(sample = c("sample 1 +", "sample 1 -",
                      "sample 2 +", "sample 2 -",
                      "sample 3 +", "sample 3 -",
                      "sample 4 +", "sample 4 -"
                      ),
           cl_sites = c(NROW(bw_1_plus), NROW(bw_1_minus),
                           NROW(bw_2_plus), NROW(bw_2_minus),
                           NROW(bw_3_plus), NROW(bw_3_minus),
                           NROW(bw_4_plus), NROW(bw_4_minus)
                           ),
           cl_events = c(sum(bw_1_plus$score), sum(bw_1_minus$score),
                         sum(bw_2_plus$score), sum(bw_2_minus$score),
                         sum(bw_3_plus$score), sum(bw_3_minus$score),
                         sum(bw_4_plus$score), sum(bw_4_minus$score)
           ))

crosslinks_overview_table

# calculate average events per location
crosslinks_overview_table$ratio <- with(crosslinks_overview_table, cl_events / cl_sites)

crosslinks_overview_table


#################
# make histograms of crosslink events
#################

# combine + and - strand of each sample to make one histogram per sample
bw_1 <- c(bw_1_plus, bw_1_minus)
bw_1

bw_2 <- c(bw_2_plus, bw_2_minus)
bw_3 <- c(bw_3_plus, bw_3_minus)
bw_4 <- c(bw_4_plus, bw_4_minus)

bw_1_df <- as.data.frame(bw_1)
bw_2_df <- as.data.frame(bw_2)
bw_3_df <- as.data.frame(bw_3)
bw_4_df <- as.data.frame(bw_4)

# make a histogram
ggplot(bw_1_df, aes(x=score))+
  geom_histogram(binwidth = 1) # with binwidth you can set the size of the bins


# Maximum number of crosslink events in this sample
max(bw_1_df$score)

# How many crosslink sites harbour just 1 crosslink event
sum(bw_1_df$score == 1)

# histogram with zoom into distribution between 1 and 20
hist_bw_1 <- ggplot(bw_1_df, aes(x=score))+
  geom_histogram(binwidth = 1)+
  coord_cartesian(xlim = c(0,20))

hist_bw_1

hist_bw_2 <- ggplot(bw_2_df, aes(x=score))+
  geom_histogram(binwidth = 1)+
  coord_cartesian(xlim = c(0,20))

hist_bw_3 <- ggplot(bw_3_df, aes(x=score))+
  geom_histogram(binwidth = 1)+
  coord_cartesian(xlim = c(0,20))

hist_bw_4 <- ggplot(bw_4_df, aes(x=score))+
  geom_histogram(binwidth = 1)+
  coord_cartesian(xlim = c(0,20))

hist_bw_1
hist_bw_2
hist_bw_3
hist_bw_4

#################
# merge bigWig files (not to be run in the course)
#################

## ## make lists of what should be combined
## # bw_all_plus <- c(bw_1_plus,bw_2_plus, bw_3_plus, bw_4_plus)
## # bw_all_minus <- c(bw_1_minus,bw_2_minus, bw_3_minus, bw_4_minus)
## #
## ## combine samples and add scores: unique_granges() from gintools package
## # bw_all_plus <- unique_granges(bw_all_plus, sum.cols = "score")
## # bw_all_minus <- unique_granges(bw_all_minus, sum.cols = "score")
## #
## ## the merged GRanges can be exported again as a bigWig file
## # path_to_save_plus <- "/path/to/your/folder/crosslinks_<RBP>_merged.plus.bw"
## # path_to_save_minus <- "/path/to/your/folder/crosslinks_<RBP>_merged.minus.bw"
## #
## # export.bw(bw_all_plus, path_to_save_plus, fixedSummaries= TRUE)
## # export.bw(bw_all_minus, path_to_save_minus, fixedSummaries= TRUE)

# import merged bigWig files
import.bw(bw_plus_path)
import.bw(bw_minus_path)

## ##########################
## # peak calling with PureCLIP (to be run outside R!)
## #######################


## #### ! do not run !!!

## pureclip -i iCLIP_reads_<RBP>_merged.bam # bam file with the merged replicates of your experiments
## -bai iCLIP_reads_<RBP>_merged.bam.bai # corresponding bai index
## -o peakcalling_<RBP>_sites.bed # name of the output sites
## -or peakcalling_<RBP>_regions.bed # a second type of output that we do not use here
## -nt 10 # nodes for calculation
## -iv 'chr1;chr2;chr3;' # chromosomes used to train the hidden markov model
##                       # (more take longer to calulate, but might be more accurate;
##                       # usually 3 chromosomes are sufficient)


#########################
# load PureCLIP sites as GRanges
#########################

pureclip_sites <- import(pureclip_path, format = "bedgraph")
pureclip_sites

# note that that strand and names of some metadata columns got lost during import
# restore
strand(pureclip_sites) <- pureclip_sites$NA.1 
mcols(pureclip_sites) <- DataFrame(score = pureclip_sites$NA.)
pureclip_sites$round_score <- round(pureclip_sites$score, digits = 1)

pureclip_sites <- keepStandardChromosomes(pureclip_sites, pruning.mode = "coarse")
pureclip_sites


#########################
# binding site definition
#########################

## Step 1: merge PureCLIP sites into regions

# merge significant crosslink sites over distances < 8 nt
pureclip = reduce(pureclip_sites, min.gapwidth = 8)

# remove sites with 1 or 2 nt length
pureclip = pureclip[width(pureclip) > 2]

# here is an example coverage vector crosslink events on subsequent positions
example_crosslinks = c(rep(0,5), rep(1,3), 0, 1, 2, rep(1,3), rep(0,4), 1, 0)
example_crosslinks
# this is the same vector in Rle encoding
Rle(example_crosslinks)

# you can obtain Rle objects for our crosslink events by applying coverage()
coverage(bw_1_plus)

## Step 2: interatively place binding sites into regions

# first we have to decide what size our binding sites should have
# as GRanges can easily be resized symetrically on both sites with +n,
# we specify the half window size
# eg. window.radius <- 3 will result in
# 3 nt downstream + 1 nt center + 3 nt upstream = a 7-nt binding site
window.size <- 3

define_binding_sites_from_matrix <- function (bw_path, pureclip, window.size){
 
  # get crosslinks as Rle   
  bw <- import.bw(bw_path, as="Rle") 
  
  # we initalise two GRanges objects 
  # 1. a GRanges objects which will contain our binding sites in the end (final BS)
  final.BS.gr <- GRanges()
  
  # 2. a GRAnges with the regions that are not yet resized (remaining regions)  
  # as in the beginning nothing is resized, we put all PureCLIP sites into remaining regions
  remaining.regions.gr <- pureclip
  
  #######################
  # loop over a matix of all PureCLIP calls to define binding sites
  #######################
  
  while(length(remaining.regions.gr) != 0){
    # get the raw CL counts in the remaining PureCLIP CL regions
    # return Rle list of all regions and turn it into matrix
    remaining.PureCLIP.CL.regions.m <- as.matrix(bw[remaining.regions.gr])
    
    # identify the center of the PureCLIP CL regions (position with max counts)
    # and store its indice
    # set NA to -Infinite
    remaining.PureCLIP.CL.regions.m[is.na(remaining.PureCLIP.CL.regions.m)] <- -Inf   
    max.pos.indice <- max.col(remaining.PureCLIP.CL.regions.m, ties.method = "first")
    
    # create a peak region that is centered to the max position
    binding_site <- remaining.regions.gr
    start(binding_site) <- start(binding_site) + max.pos.indice - 1
    width(binding_site) <- 1
    binding_site <- binding_site + window.size
    
    
    # add the new binding site in the output GRanges
    final.BS.gr <- c(final.BS.gr, binding_site)
    
    # remove the peaks from the CL regions to search for additional peaks
    # excise additionally x nucleotides up and downstream
    binding_site_surround <- as(binding_site + window.size, "GRangesList")
    
    remaining.regions.gr <- unlist(psetdiff(remaining.regions.gr, binding_site_surround))
  }
  return(final.BS.gr)
}

# define binding sites with the function above
# binding sites on plus strand
binding_sites_plus <- define_binding_sites_from_matrix(bw_path = bw_plus_path,
    pureclip = pureclip[strand(pureclip)=="+"], window.size = window.size)
# binding sites on minus strand
binding_sites_minus <- define_binding_sites_from_matrix(bw_path = bw_minus_path,
    pureclip = pureclip[strand(pureclip)=="-"], window.size = window.size)

binding_sites_plus
binding_sites_minus

# you can see that the binding sites are not sorted by chromosome anymore
# sort order of chromosomes
binding_sites_plus <- sort(binding_sites_plus)
binding_sites_minus <- sort(binding_sites_minus)
binding_sites_plus 

# merge binding sites from plus and minus strand 
Binding_sites <- c(binding_sites_plus, binding_sites_minus)

# export intermediate binding sites for visualisation in IGV
export(Binding_sites,
        con = bs_intermediate_out_file,
        format = "bed")


###########################
# keep only binding sites with at least 2 crosslink sites
############################

# make a matrix of crosslinks within binding sites
# we use an Rle of the crosslink sites
bw_plus_rle <- import.bw(bw_plus_path, as="Rle") 
bw_minus_rle <- import.bw(bw_minus_path, as="Rle") 

# each row will be one binding site
Binding_sites_crosslink_matrix_plus <- as.matrix(bw_plus_rle[binding_sites_plus])
Binding_sites_crosslink_matrix_minus <- as.matrix(bw_minus_rle[binding_sites_minus])

head(Binding_sites_crosslink_matrix_plus)

# sum up positions with crosslinks per binding site
n_crosslink_per_site_plus <- apply(Binding_sites_crosslink_matrix_plus, 
                                   1, function(x) sum(x > 0)) 
n_crosslink_per_site_minus <- apply(Binding_sites_crosslink_matrix_minus, 
                                    1, function(x) sum(x > 0)) 

head(n_crosslink_per_site_plus)

# set a cut-off
Binding_sites_n_cl_plus <- binding_sites_plus[n_crosslink_per_site_plus >= 3]
Binding_sites_n_cl_minus <- binding_sites_minus[n_crosslink_per_site_minus >= 3]

Binding_sites_n_cl <- c(Binding_sites_n_cl_plus, Binding_sites_n_cl_minus)


###########################
# add PureCLIP score of peak center
############################

# get center positions (peak summit)
bs_summits <- Binding_sites_n_cl - window.size

# overlap with PureCLIP sites
overlaps <- findOverlaps(bs_summits, pureclip_sites)
overlaps

# queryHits and subjectHits can be accessed with from() and to()
head(from(overlaps))
head(to(overlaps))

# transfer score from PureCLIP sites to binding sites
Binding_sites_n_cl$score[from(overlaps)] <- pureclip_sites$score[to(overlaps)]


###########################
# save final set of resized and filtered binding sites
############################

save(Binding_sites_n_cl, file = bs_out_file)

