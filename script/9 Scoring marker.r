library(Seurat)


# 1 Compute the 'log2FC' value after filtering by MIMOSCA and GMM
log2fc_after_filt=function(object,cell_retain,meta_label,target,nontarget){
	#@object:The subseted and proceessed Seurat object
	#@cell_retain:The list R-object of cells filtered by the function 'GMM_assignment' (8 MIMOSCA.r)
	#@meta_label:The colname indicating target(sgRNA or gene)
	#@target:The certain target(sgRNA or gene) choosed
	#@nontarget:The nontarget
	#@output:The data.frame object of results
	cell_filt=do.call(c,cell_retain[grep(names(cell_retain),pattern=target)])
	target_col=object@meta.data[,meta_label,drop=F]
	cells_target=rownames(target_col[,target_col %in% target])
	cells_target=intersect(cells_target,cell_filt)
	object_target=object[,cells_target]
	cells_nt=rownames(target_col[,target_col %in% nontarget])
	object_nt=object[,cells_nt]
	object_sub=merge(object_target,object_sub)
	object_sub=NormalizeData(object_sub)
	Idents(object_sub)=meta_label
	log2fc=FoldChange(object_sub,ident.1=target,ident.2=nontarget,mean.fxn = function(x) {
		return(log(x = rowMeans(x = expm1(x = x)) + 1, base =2))
	})
	return(log2fc)
}



# 2 Score the key markers
score_markers=function(log2fc_table,markers,target){
	#@log2fc_table:The output of the function 'log2fc_after_filt'
	#@markers:The specified markers
	#@target:The certain target(sgRNA or gene) choosed
	#@output:The scores of certain target
	log2fc_markers=log2fc_table[markers,"avg_log2FC"]
	mean_markers=round(mean(abs(log2fc_markers)),4)
	return(mean_markers)
}
 

