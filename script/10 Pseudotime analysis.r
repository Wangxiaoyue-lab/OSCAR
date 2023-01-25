library(Seurat)
library(monocle)



plot_monocle_markers=function(object,markers){
	#@object:The Seurat object
	#@markers:The markers specified
	#@output:The heatmap of pseudotime analysis
	message("Construt the CellDataSet object")
	expr_matrix <- as(as.matrix(object@assays$RNA@counts), 'sparseMatrix')
	p_data <- object@meta.data 
	f_data <- data.frame(gene_short_name = row.names(object),row.names = row.names(object))
	pd <- new('AnnotatedDataFrame', data = p_data) 
	fd <- new('AnnotatedDataFrame', data = f_data)
	cds <- newCellDataSet(expr_matrix,phenoData=pd,featureData=fd,lowerDetectionLimit=0.5,expressionFamily = negbinomial.size())
	
	message("Preprocess the CellDataSet object")
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	cds <- detectGenes(cds, min_expr = 0.1)
	expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
	object <- NormalizeData(object) 
	suppressWarnings(suppressMessages(object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE, assay = "RNA")))
	express_genes <- VariableFeatures(mouse)
	cds <- setOrderingFilter(cds, express_genes)
	plot_ordering_genes(cds)
	cds <- reduceDimension(cds, max_components = 2,method='DDRTree',num_dim=6)
	cds <- orderCells(cds)
	cds <- orderCells(cds, root_state =3)
	
	message("Draw the heatmap")
	exist_genes <- row.names(subset(fData(cds),gene_short_name %in% markers))
	Time_diff <- differentialGeneTest(cds[exist_genes,], cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
	sig_gene_names <- row.names(subset(Time_diff, qval < 0.1))
	p_heatmap=plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=1, show_rownames=T, return_heatmap=T )
	return(p_heatmap)
}

