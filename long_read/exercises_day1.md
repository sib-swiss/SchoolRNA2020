# Exercises day 1

## 1.1 First login
><img border="0" src="https://www.svgrepo.com/show/14756/person-silhouette.svg" width="30" height="30"> 30 minutes

### Login to AWS EC2 remote server
You will receive an e-mail shortly before the workshop with a key, username and IP address to login on a cloud server.
Login like this:

```sh
ssh -i path/to/key/key_[USERNAME].pem [USERNAME]@[AWS_IP]
```

### Initiate conda

Once you have logged in, initiate conda:

```sh
/opt/miniconda3/bin/conda init
```

This modifies your `.bashrc`. So logout and login again, and you can use conda.

### Setup your favourite editor to work remotely

To directly initiate and modify scripts on the remote server you can use plugins:
* Notepadd++: NppFTP
* Atom: remote-edit-ni

With the following details:
* protocol: sftp
* username: your username
* hostname: server IP
* port: 22
* authentication/logon type: path to private key file

### Some important directories

* **home:** Your home directory you can approach directly with  typing `cd`, and if you want a specific directory within home you can use tilde `~` expansion. For example:
```sh
cd ~/flair
```
* **reads:** The reads are stored in two flavours: batch separated and batch combined. This is because the authors used two different sequence runs. The different runs have different quality scores but the samples within runs is not well balanced. Therefore for quality control we use the sequence files separated in batches at:
```sh
cd /data/reads/lrrnaseq
```
A subdirectory of this directory contains the batch combined files:
```sh
cd /data/reads/lrrnaseq/batch_combined
```

* **references:** as the reference sequence we use chromosome 12 of GRCh38.p13, with GRCh38.100 as gtf.
The sequence of chr12 can be found here:
```
/data/references/GRCh38.p13.chr12.fa
```
And the gtf here:
```
/data/references/Homo_sapiens.GRCh38.100.gtf
```

## 1.2 Quality control
><img border="0" src="https://www.svgrepo.com/show/14756/person-silhouette.svg" width="30" height="30"> 1 hour

Check out the summary statistics and visualisations of a single fastq file:

```sh
NanoPlot \
--fastq /data/reads/lrrnaseq/cerebellum-5238-batch2.fastq.gz \
--N50 \
--prefix parietal_cortex-5238-batch1_ \
-o ~/nanoplot
```

### Question 1.2A:
* There is not a wide distribution of read length. Is that expected?

Also calculate the summary statistics with `fastqc` and compare them with the `NanoPlot` result:

```sh
fastqc \
/data/reads/lrrnaseq/cerebellum-5238-batch2.fastq.gz \
-o ~/fastqc
```

### Question 1.2B
* There's a pretty obvious difference in mean read quality scores between the `NanoPlot` and `fastqc` output. How does that affect the expected error rate?

> **HINT:** Find the equation to calculate error probability from quality score on [Wikipedia](https://en.wikipedia.org/wiki/Phred_quality_score)

> **SEE ALSO:** This [blog](https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/) of the author of NanoPlot, and this [thread](https://github.com/wdecoster/NanoPlot/issues/191). 

## 1.3 Read alignment
><img border="0" src="https://www.svgrepo.com/show/14756/person-silhouette.svg" width="30" height="30"> 1 hour

The sequence aligner [`minimap2`](https://github.com/lh3/minimap2) is specifically developed for (splice-aware) alignment of long reads. Checkout the helper `minimap2 --help` and/or the [github readme](https://github.com/lh3/minimap2).

### Question 1.3A:
* What would be the most logical parameter for our dataset to the option `-x`?

Introns can be quite long in mammals; up to a few hundred kb. Look up the CACNA1C gene at the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu) and roughly estimate the length of the longest intron.

### Question 1.3B:
* How does this relate to the default parameter to option `-G` of `minimap2`?

Before running `minimap2` first activate its conda environment:

```sh
conda activate flair_env
```

Then, modify the command below for `minimap2` and run it from a script.

```sh
#!/usr/bin/env bash

minimap2 \
-a \
-x [PARAMETER] \
-G [PARAMETER] \
/data/references/GRCh38.p13.chr12.fa \
/data/reads/lrrnaseq/parietal_cortex-5238-batch1.fastq.gz \
| samtools sort \
| samtools view -bh > ~/parietal_cortex-5238-batch1.bam

samtools index ~/parietal_cortex-5238-batch1.bam
```

## 1.4 Read alignment entire dataset
><img border="0" src="https://www.svgrepo.com/show/220819/group-team.svg" width="30" height="30"> 30 minutes

Start the job for the read alignment for all the fastq files. It will take a few hours, so make sure you run it in the background with `nohup` (see README.md)

```sh
#!/usr/bin/env bash

mkdir ~/read_alignment

minimap2 \
-a \
-x [PARAMETER] \
-G [PARAMETER] \
/data/references/GRCh38.p13.chr12.fa \
/data/reads/lrrnaseq/batch_combined/*.fastq.gz \
| samtools sort \
| samtools view -bh > ~/read_alignment/CACNA1C_combined.bam

samtools index ~/read_alignment/CACNA1C_combined.bam
```

### Go to:
* [Exercises day 2](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/exercises_day2.md)
* [Main page](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq)
