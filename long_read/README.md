# NCCR SIB Summer School workshop: lrRNAseq

## Pre-course preparations

We will be mainly working on an Amazon Web Services ([AWS](https://aws.amazon.com/]))  Elastic Cloud (EC2) server. Our Ubuntu server behaves like a 'normal' remote server, and can be approached through `ssh` with a username, key and IP address. All participants will we granted access to a personal home directory.

Before the course, make sure you can comfortably work on a remote server. This means that you can approach it through the shell, modify scripts and transfer files:

* SSH and scripting for windows users: [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") and [Notepad++](https://notepad-plus-plus.org/downloads/)
* SSH and scripting for mac os/linux users: [Atom](https://atom.io/) with packages like: `platform-ide-terminal` and `remote-edit-ni`
* Transferring files (both windows and unix-based systems): [FileZilla](https://filezilla-project.org/)
* Or any other preferred way to login to a remote, transfer files and work with scripts.

### Other software to install on your PC

* [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/)

### About the workshop
The general aim of the workshop is to get an introduction on:
* Assessing quality of long-read sequencing data
* Methods for sequence alignment of long reads and quantification of expression
* Differential isoform expression analysis

We will be working with data from:

Clark, M. B. et al (2020). **Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain**. Molecular Psychiatry, 25(1), 37–47. https://doi.org/10.1038/s41380-019-0583-1

The authors used full-transcript amplicon sequencing with Oxford Nanopore Technology of CACNA1C, a gene associated with psychiatric risk.

We will be performing a differential isoform expression analysis on this publicly available dataset starting from the fastq files, with the following tools:
  * Quality control: [`NanoPlot`](https://github.com/wdecoster/NanoPlot)
  * Sequence alignment: [`minimap2`](https://github.com/lh3/minimap2)
  * Isoform detection and quantification: [`flair`](https://github.com/BrooksLabUCSC/flair)
  * Differential expression analysis: [`DESeq2`](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## Excersises per day can be found here:
* [Excercises day 1](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/exercises_day1.md)
* [Excercises day 2](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/exercises_day2.md)
* [Excercises day 3](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/exercises_day3.md)
* Presentations (after the course)

## Working on the cloud server

### Some warnings

The cloud server is a temporary instance for this workshop only. Altough the computational resources should be more than enough, **it's a basic Ubuntu server, and there are no hard limits on memory or CPU usage.**
Take therefore into account that great power comes with great responsibility. Overloading it can result in a reboot, cancelling all running calculations.

### RStudio server

You can work with RStudio on the cloud server with RStudio server. You can approach it with typing in your browser: `https://[AWS_IP]:8787`. Login with your username and password. The default working directory is your home directory on the server.

### Use conda

Most software is pre-installed on the server. However, if you need specific software or versions, you can install it with `conda`. Read more on the use of `conda` [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html).

The conda base environment can not be modified, so make your own:

```sh
conda create --name my_env
```

Activate it:

```sh
conda activate my_env
```

And install whatever software you need:

```sh
conda install -c bioconda samtools=1.9
```

### Visualisation of results

Many results come in an image (e.g. `.png`, `.jpg`) or `html` format. These can not be viewed directly from the server. Therefore, transfer them to your local PC first (with e.g. FileZilla or `scp`).

### Detaching a job

On this server, there is no job scheduler, so everything is run directly from the command line. That means that if a process is running, the command line will be busy, and the job will be killed upon logout. To circumvent this, there are several methods to 'detach' the screen or prevent a 'hangup signal' of a job runnig in the background that will terminate your running job.
The software `screen` or `tmux` can be used to detach your screen, and all messages to stderr or stdout (if not redirected) will be printed to the (detached) console. Use those if you're comfortable with them.

Another, more basic, program to prevent the 'hangup signal' is `nohup`. Use it like so:

```sh
nohup [YOUR COMMAND] &
```

So, for running e.g. a shell script this would be:

```sh
nohup script.sh &
```

Anything written to stdout or stderr will be written to the file `nohup.out` in your current working directory.
