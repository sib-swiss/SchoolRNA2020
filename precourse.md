# <img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="40" height="40"> Pre-course Materials

***

<br/>

Since each of the 4 group themes have completely different softwares, please follow in instructions specific to each one:

## RNA-seq hands-on (Day 1)
***

- [Instructions to set up everything for computation](http://genoweb.toulouse.inra.fr/~sdjebali/courses/SIB_august2020/instructions/0.setup/)
- [Instructions for running the GENE-SWitCH RNA-seq pipeline](http://genoweb.toulouse.inra.fr/~sdjebali/courses/SIB_august2020/instructions/1.pipeline/)
- [Instructions for running the RNA-seq downstream analyses](http://genoweb.toulouse.inra.fr/~sdjebali/courses/SIB_august2020/instructions/2.analyses/)

<br/>

<br/>

# Projects


### <img border="0" src="/SchoolRNA2020/logos/single_cell.png" width="40" height="40"> Single cell RNA analysis
***

1. For the single cell section, we will be using Conda for managing and installing the necessary software. please follow the instructions depicted here on how to install and use conda environments: [conda instructions](conda_instructions.md). PS.: As mentioned in the instructions, Windows users are required to proceed with the Linus subsystem or via VirtualBox.

2. After installation of conda you can use one of the following files to create your environemnt for the course:
- For MacOS: [environment_macos.yml](single_cell/code/environment_macos.yml)
- For Linux/Windows: [environment_linux.yml](single_cell/code/environment_linux.yml)

The only difference between the environments is the list of compiler packages that are specific for each operationsl system. After sucessful creation of the environment following the [conda instructions](conda_instructions.md), you will be able to use 'Rstudio' with 'R 3.6.1' and 'Seurat 3.2.0' and several other associated packages.

3. Additional, but **optional** installation of other packages can be done within the environemnt via standard R packages already included in the environment ('devtools', 'install.packages' and 'biocmanager'). Some of those additional methods are: 'seurat-wrappers', 'harmony', etc.

<br/>

### <img border="0" src="/SchoolRNA2020/logos/long_read.png" width="40" height="40"> Long-read RNA sequencing
***

- [Pre-course instructions](https://github.com/GeertvanGeest/NCCR_SIB_lrRNAseq/blob/master/README.md)

<br/>

### <img border="0" src="/SchoolRNA2020/logos/ribo_profiling.png" width="40" height="40"> Ribosome-profiling
***

<br/>

### <img border="0" src="/SchoolRNA2020/logos/uv_crosslink_ip.png" width="40" height="40"> UV cross-linking immunoprecipitation (CLIP-seq)
***

<br/>

### [Back to main](README.md)
