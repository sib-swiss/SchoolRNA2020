#!/usr/bin/env bash
#
NanoPlot \
--fastq /data/reads/lrrnaseq/cerebellum-5238-batch2.fastq.gz \
--N50 \
--prefix parietal_cortex-5238-batch1_ \
-o ~/nanoplot

fastqc \
/data/reads/lrrnaseq/cerebellum-5238-batch2.fastq.gz \
-o ~/fastqc

## approx. 10 minutes on 3 cores
minimap2 \
-a \
-x splice \
-G 500k \
/data/references/GRCh38.p13.chr12.fa \
/data/reads/lrrnaseq/parietal_cortex-5238-batch1.fastq.gz \
| samtools sort \
| samtools view -bh > ~/parietal_cortex-5238-batch1.bam

samtools index ~/parietal_cortex-5238-batch1.bam

## start: 11:11
mkdir ~/read_alignment

## takes about three hours..
minimap2 \
-a \
-x splice \
-G 500k \
/data/references/GRCh38.p13.chr12.fa \
/data/reads/lrrnaseq/batch_combined/*.fastq.gz \
| samtools sort \
| samtools view -bh > ~/read_alignment/CACNA1C_combined.bam

samtools index ~/read_alignment/CACNA1C_combined.bam

# mkdir ~/flair_output
#
python3 ~/flair/bin/bam2Bed12.py \
-i ~/read_alignment/CACNA1C_combined.bam \
> ~/flair_output/CACNA1C_combined.bed12
#
python3 ~/flair/flair.py correct \
-q ~/flair_output/CACNA1C_combined.bed12 \
-g /data/references/GRCh38.p13.chr12.fa \
-f /data/references/Homo_sapiens.GRCh38.100.gtf \
-o ~/flair_output/CACNA1C
#
READS=`ls /data/reads/lrrnaseq/batch_combined/*.fastq.gz | tr "\n" ","`
READS="${READS%?}" #remove last comma

python3 ~/flair/flair.py collapse \
-g /data/references/GRCh38.p13.chr12.fa \
-r $READS \
-q ~/flair_output/CACNA1C_all_corrected.bed \
-f /data/references/Homo_sapiens.GRCh38.100.gtf \
-s 3 \
-t 3 \
-o ~/flair_output/CACNA1C.collapse
#
cd ~/flair_output

python3 ~/flair/flair.py quantify \
-r /data/reads/lrrnaseq/batch_combined/reads_manifest.tsv \
-i ~/flair_output/CACNA1C.collapse.isoforms.fa \
-t 4
