# <img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="40" height="40"> Pre-course Materials

***

<br/>

Here is a list of instructions for installing packages and other required software for the course.

- [Install Docker](docker_instructions.md)

Since each of the 4 group themes have completely different softwares, please follow in instructions specific to each one:

## RNA-seq hands-on

For the RNA-seq hands-on there are two possibilities :
- either you install singularity=3.5.3 (see instructions here: https://sylabs.io/guides/3.5/user-guide/quick_start.html, be careful you need to be root to do it) 
  - and then you use the singularity image we provide (you may have to install some other tools but it will be quick)

- or you install the following tools and version :
  - python=3.7.4
  - R=3.5.1 (with the following libraries: optparse, reshape2, ggplot2, plyr, GOstats, EdgeR)
  - nextflow=20.01.0  (instructions here: https://www.nextflow.io/docs/latest/getstarted.html)
  - fastqc=0.11.7
  - trimgalore=0.6.5
  - cutadapt=2.10
  - fastqc=0.11.9
  - pandas=1.0.3
  - samtools=1.10
  - star=2.7.3a
  - stringtie=2.1.1


## Projects


### single cell RNA analysis

- [Install conda](conda_instructions.md)


### long-read RNA sequencing


### ribosome-profiling


### UV cross-linking immunoprecipitation (CLIP-seq)



<br/>

### [Back to main](README.md)
