# Exercises day 2

## 2.1 Prepare FLAIR
><img border="0" src="https://www.svgrepo.com/show/14756/person-silhouette.svg" width="30" height="30"> 30 minutes

FLAIR is a set of python scripts that can be used to identify and quantify (new) isoforms based on alignment files of long-read sequencing data. We're basically following the pipeline as described [here](https://github.com/BrooksLabUCSC/flair).

To use FLAIR first clone the git repository, and activate the (pre-configured) conda environment:

```sh
cd
git clone https://github.com/BrooksLabUCSC/flair.git
conda activate flair_env
```

After that, generate a FLAIR working directory and convert our `.bam` alignment file to `bed12` format.

```sh
mkdir ~/flair_output

python3 ~/flair/bin/bam2Bed12.py \
-i ~/read_alignment/CACNA1C_combined.bam \
> ~/flair_output/CACNA1C_combined.bed12
```

## 2.2 Run flair
><img border="0" src="https://www.svgrepo.com/show/220819/group-team.svg" width="30" height="30"> 2 hours

Generate (and run) a shell script to run the modules `flair.py correct`, `flair.py collapse` and `flair.py quantify`. To do this, carefully follow the manual at https://github.com/BrooksLabUCSC/flair. Structure and document your script(s), so you can easily re-run the analysis.

Files you will need are:
* Reference genome (chromosome 12 only): `/data/references/GRCh38.p13.chr12.fa`
* GTF: `/data/references/Homo_sapiens.GRCh38.100.gtf`
* Reads manifest: `/data/reads/lrrnaseq/batch_combined/reads_manifest.tsv`

> **HINT1:** To get all fastq files comma-separated in a variable (for `flair.py collapse`), you can use this code:

```sh
READS=`ls /data/reads/lrrnaseq/batch_combined/*.fastq.gz | tr "\n" ","`
READS="${READS%?}" #remove last comma
```
> **HINT2:** For visualisation of (intermediate) files you can use IGV. It accepts `.bam`, `.bed12` and `.gtf` files. The human genome hg38 should already be available after you installed IGV.

### Go to:
* [Exercises day 3](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/exercises_day3.md)
* [Main page](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq)
