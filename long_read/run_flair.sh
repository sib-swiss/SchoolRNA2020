#!/bin/bash

# mkdir ~/flair_output
# mkdir ~/logs

## minimap2 runs in about 90 mins
## probably less time with 6 cpu

minimap2 \
-ax splice \
-G 500k \
/data/lrrnaseq/references/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
/data/lrrnaseq/reads/batch_combined/*.fastq.gz \
2> ~/logs/minimap2.log \
| samtools sort \
| samtools view -bh > ~/flair_output/CACNA1C_combined.bam

samtools index ~/flair_output/CACNA1C_combined.bam

# use hg38 in IGV for visualisation. works.

#
python3 ~/flair/bin/bam2Bed12.py \
-i ~/flair_output/CACNA1C_combined.bam \
> ~/flair_output/CACNA1C_combined.bed12
#
python3 ~/flair/flair.py correct \
-q ~/flair_output/CACNA1C_combined.bed12 \
-g /data/reads/lrrnaseq/references/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-f /data/reads/lrrnaseq/references/Homo_sapiens.GRCh38.100.gtf \
-o ~/flair_output/CACNA1C 2> ~/logs/flair_correct.log


READS=`ls /data/reads/lrrnaseq/reads/batch_combined/*.fastq.gz | tr "\n" ","`
READS="${READS%?}" #remove last comma

python3 ~/flair/flair.py collapse \
-g /data/reads/lrrnaseq/references/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-r $READS \
-q ~/flair_output/CACNA1C_all_corrected.bed \
-f /data/reads/lrrnaseq/references/Homo_sapiens.GRCh38.100.gtf \
-s 3 \
-t 6 \
-o ~/flair_output/CACNA1C.collapse \
--keep_intermediate \
2> ~/logs/flair_collapse.log

cd ~/flair_output

python3 ~/flair/flair.py quantify \
-r ~/reads_manifest.tsv \
-i ~/flair_output/CACNA1C.collapse.isoforms.fa \
-t 6 \
2> ~/logs/flair_quantify.log

python3 ~/flair/flair.py diffExp \
-q ~/flair_output/counts_matrix.tsv \
-o ~/flair_output/diffexp \
-t 4 \
2> ~/logs/flair_diffExp.log


# python3 ./flair/flair.py diffSplice \
# -i CACNA1C.collapse.isoforms.bed \
# -q counts_matrix.tsv \
# 2> flair_diffSplice.log

python3 ./flair/bin/predictProductivity.py \
-i CACNA1C.collapse.isoforms.bed \
-g ./reference/Homo_sapiens.GRCh38.100.gtf \
-f ./reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
--longestORF \
2> predictProductivity.log \
> CACNA1C_productivity.bed

# export PATH=$PATH:/Users/geertvangeest/Documents/NCCR_summerschool/flair:/Users/geertvangeest/Documents/NCCR_summerschool/flair/bin
# cd ./flair/bin/
#
# python3 plot_isoform_usage.py \
# ../../CACNA1C_productivity.bed \
# ../../counts_matrix.tsv \
# ENSG00000151067
