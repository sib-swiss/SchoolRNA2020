




library(future)
plan(multiprocess, workers=(availableCores()-1) )

PATH <- "~/Downloads/Marrow/"

#######
# Parallelization for reading input files from a folder and creating a count matrix
#######
files <- list.files(PATH)
genes <- rownames(read.csv(paste0(PATH,files[1]),header = F,row.names = 1))


file_chunk <- cut(as.numeric(factor(files)),breaks = 7)
file_chunk <- lapply( levels(file_chunk) , function(x){ files[file_chunk==x] })

f <- list()
for(ind in 1:length(file_chunk) ){
  f[[ ind ]] <- future({
    rrr <- lapply(file_chunk[[ind]],function(x){
      return(Matrix::Matrix(read.csv(paste0(PATH,x),header = F,row.names = 1)[,1],sparse = T))})
    rrr <- do.call(cbind,rrr)
    colnames(rrr) <- sub("[.]merge.*","",file_chunk[[ind]])
    return(rrr)
  })
}
#######




#######
# Compiles the results from each core into one single matrix
#######
metadata_all <- lapply(f, FUN = value)
metadata_all <- do.call(cbind,metadata_all)
dim(metadata_all)
rownames(metadata_all) <- genes
# Matrix::writeMM(metadata_all,file = "~/Downloads/Marrow/Marrow.mtx")
#######



library(rhdf5)
# data <- read.csv("~/Downloads/Marrow.h5",row.names = 1)
# data <- Matrix::Matrix(as.matrix(data), sparse = T)
data <- metadata_all

h5createFile("~/Downloads/Marrow.h5")
h5createGroup("~/Downloads/Marrow.h5","matrix")

h5write(data@Dimnames[[2]],"~/Downloads/Marrow.h5","matrix/barcodes")
h5write(data@x,"~/Downloads/Marrow.h5","matrix/data")
h5write(data@i,"~/Downloads/Marrow.h5","matrix/indices")
h5write(data@p,"~/Downloads/Marrow.h5","matrix/indptr")
h5write(data@Dim,"~/Downloads/Marrow.h5","matrix/shape")

h5createGroup("~/Downloads/Marrow.h5","matrix/features")
h5write(data@Dimnames[[1]]
,"~/Downloads/Marrow.h5","matrix/features/name")
h5write(data@Dimnames[[1]]
        ,"~/Downloads/Marrow.h5","matrix/features/_all_tag_keys")
h5write(rep("expression",nrow(data))
        ,"~/Downloads/Marrow.h5","matrix/features/feature_type")
h5write(rep("mm10",nrow(data))
        ,"~/Downloads/Marrow.h5","matrix/features/genome")


h5ls("~/Downloads/Marrow.h5")

nd <- Seurat::Read10X_h5("~/Downloads/Marrow.h5")
sum(!nd == data)





