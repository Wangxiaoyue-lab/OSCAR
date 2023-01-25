rm(list=ls())
library(Seurat)
getwd()


# 1 The Seurat object (proceessed) of pilot experiment
##@input:The file of expression (UMI counts) and meta.data
##@output:The Seurat object
exp=read.csv("./pilot_expression_matrix.csv.gz",header=T,row.names=1)
meta.data=read.csv("./pilot_metadata.csv",header=T,row.names=1)
rownames(meta.data)=gsub("-",".",rownames(meta.data))
identical(colnames(exp),rownames(meta.data))
DM_EM <- CreateSeuratObject(counts = exp)
DM_EM <- AddMetaData(object = DM_EM, metadata = meta.data)
save(DM_EM,file="./pilot_DM_EM.Rdata")
#load(file="./pilot_DM_EM.Rdata")



# 2 The Seurat object (proceessed) of OSCAR
##@input:The file of expression (UMI counts) and meta.data
##@output:The Seurat object
DM=read.csv("./OSCAR_DM_expression_matrix.csv.gz",header=T,row.names=1)
EM=read.csv("./OSCAR_EM_expression_matrix.csv.gz",header=T,row.names=1)
exp=cbind(DM,EM)
meta.data=read.csv("./OSCAR_metadata.csv",header=T,row.names=1)
rownames(meta.data)=gsub("-",".",rownames(meta.data))
identical(colnames(exp),rownames(meta.data))
data.filt <- CreateSeuratObject(counts = exp)
data.filt <- AddMetaData(object = data.filt, metadata = meta.data)
save(data.filt,file="./filt_dm_em2.Rdata")
#load(file="./filt_dm_em2.Rdata")


# 3 The matrix for SCENIC
write.csv(t(as.matrix(data.filt@assays$RNA@counts)),file = "./data.filt.sample.csv")
