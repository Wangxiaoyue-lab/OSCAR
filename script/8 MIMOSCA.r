library(Seurat)
library(nnet)
library(mixtools)
source("https://github.com/quadbiolab/organoid_regulomes/blob/main/crop_seq/ko_inference.R")



# 1 Extract the Y matrix
#@input:The Seurat object
load(file="./filt_dm_em2.Rdata")
data_sub = subset(x=data.filt,subset=type=="DM")
data_sub = NormalizeData(data_sub)
data_sub = CellCycleScoring(object = data_sub, g2m.features = stringr::str_to_title(cc.genes$g2m.genes),s.features = stringr::str_to_title(cc.genes$s.genes))
suppressWarnings(suppressMessages(data_sub = FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE, assay = "RNA")))



# 2 Extract and embed the X matrix
#@input:The meta.data of Seurat object
#@output:The Seurat object added the data of sgRNA
meta=data_sub@meta.data
meta$barcode=gsub("non_target","nontarget",meta$barcode)
sgrna=class.ind(meta$barcode)
rownames(sgrna)=rownames(meta)
data_sub = add_guide_assay.Seurat(object=data_sub,guide_assignments=sgrna,assay_name = 'guide_assignments')



# 3 Infer the probability of being perturbed
data_sub = infer_ko_probs.Seurat(object=data_sub,covariates = c("nFeature_RNA","Phase"),guide_sep = '-',assay = 'RNA',slot = 'data',alpha = 0.5,assay_name = 'perturb',guide_assay = 'guide_assignments',save_coefs = TRUE,genes_use = NULL,verbose = TRUE
)



# 4 GMM
GMM_assignment=function(guide_assignments,perturb,sgRNAs,p_threshold){
	#@guide_assignments:The 二分类化 matrix of guide assignment
	#@perturb:The matrix of probabilities of perturbations
	#@sgRNAs:The selected sgRNAs (i.e. sgRNAs of P2 and P4)
	#@p_threshold:The threshold of probabilities for cells retaining 
	#@output:
	assign_sgRNAs=guide_assignments[sgRNAs,]
	perturb_sgRNAs=perturb[sgRNAs,]
	cell_retain=list()
	for(i in 1:nrow(assign_sgRNAs){
		assign_i=assign_sgRNAs[i,]
		sgrna_i=rownames(assign_sgRNAs)[i]
		prob_i=perturb_sgRNAs[i,]
		if(sum(assign_i)>40){
			prob_i_1=prob_i[which(assign_i>0)]
			result_EM = mixtools::normalmixEM(prob_i_1, k=2)
			posterior=as.data.frame(result_EM$posterior)
			delta_max=prob_i_1[which(posterior$comp.2==max(posterior$comp.2))]-prob_i_1[which(posterior$comp.1==max(posterior$comp.1))]
			if(minus_max>0){posterior=posterior[,c(2,1)]}else{posterior=posterior[,c(1,2)]}
			colnames(posterior)=c("high","low")
			cell_high=prob_i_1[which(posterior$high > p_threshold)]
			cell_low=prob_i_1[which(posterior$high <= p_threshold)]
			delta_mean=mean(cell_high)-mean(cell_low)
			if(is.na(delta_mean)){delta_mean=0}
			if(as.numeric(delta_mean)>0.3){
				cell_success=names(cell_high)
			}else{
				cell_success=names(prob_i_1[prob_i_1>0.2])
			}
			cell_perturb=list(cell_success)
			names(cell_perturb) = sgrna
			cell_retain=c(cell_retain,cell_perturb)
		}else{
			message(paste0("The number of cells of",sgrna_i," is less than 40")) 
		}
	}
	return(cell_retain)
}

guide_assignments=as.matrix(data.DM$guide_assignments@data)
perturb=as.matrix(data.DM[["perturb"]]) 
cell_retain=GMM_assignment(guide_assignments=guide_assignments,perturb=perturb,sgRNAs=c("Arid3a-1","Arid3a-2","Arid3a-3"),p_threshold=0.05)