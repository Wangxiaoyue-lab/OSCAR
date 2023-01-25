library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(ggplotify)
library(pheatmap)
library(cowplot)
library(mclust)
library(reshape2)

##@input:The proceessed Seurat object of OSCAR
##@input:The result of SCENIC
load(file="./filt_dm_em2.Rdata")
loom = open_loom("./data.filt.out_SCENIC.loom")


#1 Extract the AUC matrix
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons") 
regulons = regulonsToGeneLists(regulons_incidMat) 
regulonAUC = get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
AUC_value=regulonAUC@assays@data@listData$AUC
meta_OSCAR=data.filt@meta.data[data.filt@meta.data$type=="DM",]
count_sgrna=as.data.frame(table(meta_OSCAR$barcode))
count_sgrna=count_sgrna[count_sgrna$Freq>40,]$Var1
meta_OSCAR=meta_OSCAR[meta_OSCAR$barcode %in% count_sgrna,]



#2 Construt the regulons X Perturbations matrix
mean_target=function(meta,AUC,target){
	#@output:The regulons X Perturbations matrix
	#@meta:The meta.data of Seurat object
	#@AUC:The matrix of AUC value
	#@target:The choice to calulate the mean(gene or sgRNA) 
	AUC=as.data.frame(AUC[,match(rownames(meta),colnames(AUC))]) 
	colnames(meta)[colnames(meta) %in% target]="group"
	count_sgrna=as.data.frame(table(meta$group))
	count_sgrna=count_sgrna[count_sgrna$Freq>=40,]$Var1 
	for(i in 1:length(count_sgrna)){
		sgrna_i=count_sgrna[i]
		cell_names=rownames(meta[meta$group %in% sgrna_i,])
		AUC_i=AUC[,colnames(AUC) %in% cell_names]
		AUC_i = apply(AUC_i,1:2,function(x){as.numeric(x)})
		for(j in 1:nrow(AUC_i)){
			#value_j=shapiro.test(as.numeric(AUC_i[j,]))$p.value
			value_j=mean(as.numeric(AUC_i[j,]))
			if(j==1){AUC_i_value=value_j}else{AUC_i_value=c(AUC_i_value,value_j)}
		}
		if(i==1){mean_regulon=AUC_i_value}else{mean_regulon=cbind(mean_regulon,AUC_i_value)}
	}
	rownames(mean_regulon)=rownames(AUC)
	colnames(mean_regulon)=count_sgrna
	return(mean_regulon)
}
mean_target_OSCAR=mean_target(meta=meta_OSCAR,AUC=AUC_value,target="gene") 



#3 Obtain the effect size (perturbation - NTC)
non_target=mean_target_OSCAR[,c(45)] #Non-targeting
effect_matrix=mean_target_OSCAR[,-c(45)]
for(i in 1:nrow(effect_matrix)){
	non_target_i=non_target[i]
	effect_matrix_i=effect_matrix[i,]-non_target_i
	effect_matrix[i,]=effect_matrix_i
}




#4 Draw the correlation heatmap and cluster the perturbations/regulons
kmeans_heatmap=function(matrix,per_k,reg_k){
	#@matrix:The matrix of effect size
	#@per_k:The parameter k for perturbation
	#@reg_k:The parameter k for regulon
	#@output:p_grid--The combined heatmap
	#@output:corper_k_gene--The result of kmeans for perturbation
	#@output:correg_k_gene--The result of kmeans for regulon
	#4.1 the heatmap before rearrangement
	raw=t(scale(t(matrix)))
	raw[raw > 4] = 4 
	raw[raw < -4]= -4
	p_raw=as.ggplot(pheatmap(raw,show_colnames =F,show_rownames = F, border = F, color = colorRampPalette(colors = c("navy","white","red"))(100)))

	#4.2 the correlation heatmap of perturbation
	cor_per=as.data.frame(cor(raw))
	cor_per[cor_per > 0.5] = 0.5 
	cor_per[cor_per < -0.5]= -0.5
	cor_per_kmeans=kmeans(cor_per,per_k)
	corper_k=sort(cor_per_kmeans$cluster)
	corper_k_gene=as.data.frame(corper_k)
	colnames(corper_k_gene)[1] = "perturb_kmeans"
	corper_k_gene$gene=rownames(corper_k_gene)
	corper_arrange=names(corper_k)
	cor_per=cor_per[match(corper_arrange,colnames(cor_per)),match(corper_arrange,rownames(cor_per))]
	group_color = c("#85B22E","#5F80B4","#E29827","#922927",'#57C3F3') 
	names(group_color) = c(1,2,3,4,5)
	anno_colors = list(perturb_kmeans=group_color)
	p_per=as.ggplot(pheatmap(cor_per,show_colnames =F,show_rownames = F, border = F, color = colorRampPalette(colors = c("purple","white","green"))(100),cluster_rows=FALSE,cluster_cols=FALSE,annotation_col = corper_k_gene,annotation_legend=FALSE,legend_labels = FALSE,annotation_names_row = FALSE, annotation_names_col = FALSE,annotation_colors = anno_colors))

	#4.3 the correlation heatmap of regulon
	cor_reg=as.data.frame(cor(t(raw)))
	cor_reg[cor_reg > 0.5] = 0.5 
	cor_reg[cor_reg < -0.5]= -0.5
	cor_reg_kmeans=kmeans(cor_reg,reg_k)
	correg_k=sort(cor_reg_kmeans$cluster)
	correg_k_gene=as.data.frame(correg_k)
	colnames(correg_k_gene)[1] = "regulon_kmeans"
	correg_k_gene$gene=rownames(correg_k_gene)
	correg_arrange=names(correg_k)
	cor_reg=cor_reg[match(correg_arrange,colnames(cor_reg)),match(correg_arrange,rownames(cor_reg))]
	group_color = c("#708090",'#68A180','#F3B1A0', '#D6E7A3','#57C3F3') 
	names(group_color) = c(1,2,3,4,5)
	anno_colors = list(regulon_kmeans=group_color)
	p_reg=as.ggplot(pheatmap(cor_reg,show_colnames =F,show_rownames = F, border = F, color = colorRampPalette(colors = c("purple","white","green"))(100),cluster_rows=FALSE,cluster_cols=FALSE, annotation_row = correg_k_gene,annotation_legend=FALSE,legend_labels = FALSE,annotation_names_row = FALSE, annotation_names_col = FALSE,annotation_colors = anno_colors))

	#4.4 the heatmap after rearrangement
	rearrange = raw[match(correg_arrange,rownames(matrix)),match(corper_arrange,colnames(matrix))]
	rearrange[rearrange > 3] = 3
	rearrange[rearrange < -3]= -3
	p_after=as.ggplot(pheatmap(rearrange,show_colnames = F,show_rownames = F, border = F, color = colorRampPalette(colors = c("navy","white","red"))(100),cluster_rows=FALSE,cluster_cols=FALSE))
	p_grid=plot_grid(p_raw,p_per, p_reg,p_after, ncol = 2)
	return(list(plot=p_grid, kmeans_perturbation=corper_k_gene,kmeans_regulon=correg_k_gene))
}

result_kmeans=kmeans_heatmap(matrix=effect_matrix,per_k=4,reg_k=5)
p_heatmap=result_kmeans[["plot"]]
kmeans_perturbation=result_kmeans[["kmeans_perturbation"]]
kmeans_regulon=result_kmeans[["kmeans_regulon"]]




#5 wilcox rank sum test
wilcox_filt=function(matrix,metadata,nt_sgrnas,AUC,Level,target){
	#@matrix:The matrix of effect size
	#@metadata:The annotation of cells (i.e. metadata from Seurat)
	#@nt_sgrnas:The regular expression pattern of Non-targeting sgRNAs
	#@AUC:The matrix of AUC value generated from SCENIC
	#@Level:The level of wilcox rank sum test(i.e. sgRNAs or cells)
	#@target:The choice to calulate the mean(gene or sgRNA)
	#@output:The matrix of P value computed by wilcox
	AUC=as.data.frame(AUC)
	mean_sgrna_OSCAR=mean_target(meta=metadata,AUC=AUC,target=target) 
	mean_sgrna_nt=mean_sgrna_OSCAR[,grep(colnames(mean_sgrna_OSCAR),pattern=nt_sgrnas)]
	colnames(metadata)[colnames(metadata) %in% target]="group"
	pvalue_matrix=as.data.frame(matrix("",nrow=nrow(matrix),ncol=ncol(matrix)))
	colnames(pvalue_matrix)=colnames(matrix)
	rownames(pvalue_matrix)=rownames(matrix)
	if(Level=="sgRNAs"){
		for(i in 1:ncol(matrix)){
			for(j in 1:nrow(matrix)){
				p=wilcox.test(as.numeric(mean_sgrna_nt[j,]), mu = as.numeric(matrix[j,i]), alternative = "two.sided")$`p.value`
				pvalue_matrix[j,i]=round(p,5)
			}
		}
	}else if(Level=="cells"){
		for(i in 1:ncol(matrix)){
			target_i=colnames(matrix)[i]
			cells_i=rownames(metadata[grep(metadata$group,pattern=paste0(target_i,"-")),])
			AUC_i=AUC[,colnames(AUC) %in% cells_i]
			sample_i=sample(ncol(AUC_i),ifelse(ncol(AUC_i)<1500,ncol(AUC_i),1500))
			AUC_i=AUC_i[,sample_i]
			AUC_i=AUC_i+0.00000001
			cells_nt=rownames(metadata[grep(metadata$group,pattern=paste0(nt_sgrnas,"-")),])
			AUC_nt=auc[,colnames(AUC) %in% cells_nt]
			sample_nt=sample(ncol(AUC_i),ifelse(ncol(AUC_i)<1500,ncol(AUC_i),1500))
			AUC_nt=AUC_nt[,sample_nt]
			AUC_nt=AUC_nt+0.00000001
			for(j in 1:nrow(AUC)){
				p=wilcox.test(as.numeric(AUC_i[j,]),as.numeric(AUC_nt[j,]), alternative = "two.sided")$`p.value`
				pvalue_matrix[j,i]=round(p,5)
			}
		}
	}
	return(pvalue_matrix)
}

p_matrix=wilcox_filt(matrix=effect_matrix,metadata=meta_OSCAR,nt_sgrnas="non_target",AUC=AUC_value,Level="sgRNAs",target="barcode")




#6 hypergeometric test
hypergeometric_test=function(effect_matrix,pvalue_matrix,kmeans_perturbation,kmeans_regulon){
	#@effect_matrix:The matrix of effect size
	#@pvalue_matrix:The matrix of P value computed by wilcox
	#@kmeans_perturbation:The result of kmeans for perturbation
	#@kmeans_regulon:The result of kmeans for regulon
	#@output:The result of hypergeometric test
	effect_to_reshape = as.data.frame(t(effect_matrix))
	effect_to_reshape$perturb=rownames(effect_to_reshape)
	effect=reshape2::melt(effect_to_reshape,id.vars=c("perturb"),variable.name = "regulon")
	colnames(effect)[3]="effect"
	pvalue_to_reshape = as.data.frame(t(pvalue_matrix))
	pvalue_to_reshape$perturb=rownames(pvalue_to_reshape)
	pvalue=reshape2::melt(pvalue_to_reshape,id.vars=c("perturb"),variable.name = "regulon")
	colnames(pvalue)[3]="pvalue"
	effect=cbind(effect,pvalue[,3,drop=F])
	effect$perturb_group=apply(as.data.frame(effect$perturb),1,function(x){kmeans_perturbation[kmeans_perturbation$gene==x,]$perturb_kmeans})
	effect$regulon_group=apply(as.data.frame(effect$regulon),1,function(x){kmeans_regulon[kmeans_regulon$gene==x,]$regulon_kmeans})
	effect$direction=ifelse(effect$effect>0,"P","N")
	Regulation=effect[effect$pvalue<0.05,]
	Regulation_stastic = c()
	Regulation_P = Regulation[Regulation$direction == "P", ]
	Regulation_N = Regulation[Regulation$direction == "N", ]
	for(i in 1:length(unique(kmeans_perturbation$perturb_kmeans))) {
		Regulation_from.i = Regulation[Regulation$perturb_group == i, ]
		Regulation_from.i_P = Regulation_from.i[Regulation_from.i$direction == "P", ]
		Regulation_from.i_N = Regulation_from.i[Regulation_from.i$direction == "N", ]
		for(j in 1:length(unique(kmeans_regulon$regulon_kmeans))){
			Regulation_to.j = Regulation[Regulation$regulon_group == j, ]
			Regulation_to.j_P = Regulation_to.j[Regulation_to.j$direction == "P", ]
			Regulation_to.j_N = Regulation_to.j[Regulation_to.j$direction == "N", ]
			Regulation_from.i_to.j_P = Regulation_from.i_P[Regulation_from.i_P$regulon_group == j, ]
			Regulation_from.i_to.j_N = Regulation_from.i_N[Regulation_from.i_N$regulon_group == j, ]
			Regulation_from.i_to.j_P_num = c(i, j, "P",round(mean(as.numeric(Regulation_from.i_to.j_P$effect)),3), nrow(Regulation_from.i_to.j_P),nrow(Regulation_to.j_P), nrow(Regulation_P) - nrow(Regulation_to.j_P), nrow(Regulation_from.i_P))
			Regulation_from.i_to.j_N_num = c(i, j, "N",round(mean(as.numeric(Regulation_from.i_to.j_N$effect)),3), nrow(Regulation_from.i_to.j_N),nrow(Regulation_to.j_N), nrow(Regulation_N) - nrow(Regulation_to.j_N), nrow(Regulation_from.i_N))
			Regulation_stastic = rbind(Regulation_stastic, Regulation_from.i_to.j_P_num, Regulation_from.i_to.j_N_num)
		}
	}
	Regulation_stastic = as.data.frame(Regulation_stastic)
	Regulation_stastic = Regulation_stastic[Regulation_stastic$V4 != "NaN", ]
	Regulation_phyper = apply(Regulation_stastic[, 5:ncol(Regulation_stastic)], 1, function(x){
		x = as.numeric(x)
		if(x[1]<4){y = 1}else{y = phyper(x[1],x[2],x[3],x[4],lower.tail = FALSE)}
	})
	Regulation_padj = p.adjust(Regulation_phyper, method = "BH")
	Regulation_stastic = cbind(Regulation_stastic[,1:4],round(Regulation_phyper,3),round(-log(as.numeric(Regulation_padj),10),3))
	colnames(Regulation_stastic)[1:6] = c("perturb_module", "regulon_module", "direction","mean_effect","Pvalue", "BH")
	rownames(Regulation_stastic) = paste0(Regulation_stastic$perturb_module,"→",Regulation_stastic$regulon_module,Regulation_stastic$direction)
	return(Regulation_stastic)
}

result_hyper=hypergeometric_test(effect_matrix=effect_matrix,pvalue_matrix=p_matrix,kmeans_perturbation=kmeans_perturbation,kmeans_regulon=kmeans_regulon)
result_hyper=result_hyper[result_hyper$Pvalue<0.05,]