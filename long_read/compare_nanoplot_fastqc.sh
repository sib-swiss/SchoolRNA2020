#!/bin/bash

NanoPlot \
--fastq /data/lrrnaseq/reads/cerebellum-5238-batch2.fastq.gz \
--N50 \
--prefix cerebellum-5238-batch2_ \
-o ~/compare_nanoplot_fastqc

/home/ubuntu/bin/FastQC/fastqc \
/data/lrrnaseq/reads/cerebellum-5238-batch2.fastq.gz \
-o ~/compare_nanoplot_fastqc
