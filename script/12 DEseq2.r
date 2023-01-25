library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)

#@exp:The expressopm matrix generated from FeatureCount
#@re_order:The rearrangement of colnames
#@condition:The specified condition
#@contrast:The specified contrast

# 1 Construt the DESeqDataSet object
mat=exp[,colnames(exp) %in% reorder]
mat=exp[,match(re_order,colnames(exp))]
metadata <- data.frame(sample_id = reorder)
metadata$condition <- factor(condition, levels = unique(condition))
dds <-DESeqDataSetFromMatrix(countData=mat,colData=metadata,design=~condition,tidy=F)
dds <- dds[rowSums(counts(dds))>1,]
dim(dds)


# 2 exploratory analysis
vsd <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
hc <- hclust(sampleDists, method = "ward.D")
plot(hc, hang = -1)
plotPCA(vsd,"condition")


# 3 Differential gene analysis
dds <- DESeq(dds)
res <- results(dds,contrast=contrast, alpha = 0.1)
plotMA(dd1, ylim=c(-10,10))
resLFC <- lfcShrink(dds, contrast=contrast, type="ashr")
plotMA(resLFC, ylim=c(-10,10))
res$gene_id=rownames(res)
resLFC$gene_id=rownames(gene_id)

# 4 Transform gene ID
res$symbol <- mapIds(org.Mm.eg.db,keys=res$gene_id,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,keys=res$gene_id,column="ENTREZID",keytype="ENSEMBL",multiVals="first")

resLFC$symbol <- mapIds(org.Mm.eg.db,keys=resLFC$gene_id,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
resLFC$entrez <- mapIds(org.Mm.eg.db,keys=resLFC$gene_id,column="ENTREZID",keytype="ENSEMBL",multiVals="first")


