#!/bin/bash

mkdir ~/quality_control
mkdir ~/quality_plots

BATCH1=`ls /data/lrrnaseq/reads/*batch1.fastq.gz`
BATCH2=`ls /data/lrrnaseq/reads/*batch2.fastq.gz`

NanoPlot --fastq $BATCH1 \
--N50 \
--prefix batch1_ \
-o ~/quality_plots \
--readtype 2D

NanoPlot --fastq $BATCH2 \
--N50 \
--prefix batch2_ \
-o ~/quality_plots \
--readtype 2D

# cd /data/lrrnaseq/reads
# READF=`ls ERR*.fastq.gz`
#
# for readfile in $READF
# do
#   readID=`echo "$readfile" | cut -f 1 -d '.' `
#   mkdir ~/quality_control/$readID
#   nanoQC -o ~/quality_control/$readID \
#   /data/lrrnaseq/reads/$readfile
# done

## cd ~/quality_control/ERR3577080
## python3 -m http.server 8000
## view like (provided port 8000 is open):
## http://[SERVERIP]/nanoQC.html:8000

/home/ubuntu/bin/FastQC/fastqc \
/data/lrrnaseq/reads/cerebellum-5238-batch2.fastq.gz \
-o ~/fastqc_out

NanoPlot --fastq /data/lrrnaseq/reads/cerebellum-5238-batch2.fastq.gz \
--N50 \
--prefix cerebellum-5238-batch2_ \
-o ~/quality_plots \
--readtype 2D

NanoPlot --fastq /data/lrrnaseq/reads/cerebellum-5238-batch2.fastq.gz \
--N50 \
--prefix cerebellum-5238-batch2_1D2_ \
-o ~/quality_plots \
--readtype 1D2
