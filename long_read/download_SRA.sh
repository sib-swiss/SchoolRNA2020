#!/bin/bash


cd /data/lrrnaseq/reads

## last two entries are raw data
for SRR in `tail -n +2 SraRunInfo.csv | head -24 | cut -f 1 -d ',' `
do
  echo $SRR

  prefetch $SRR && \
  fastq-dump \
  --gzip \
  $SRR
done

## give new filenames
while read OLD NEW
do
  mv $OLD $NEW
done < ~/NCCR_SIB_lrRNAseq/file_rename.txt

## combine files from two batches
BATCH=`ls *.fastq.gz | cut -f 1,2 -d "-" | uniq`

mkdir batch_combined

for PREFIX in $BATCH
do
  cat $PREFIX* > ./batch_combined/$PREFIX.fastq.gz
done
