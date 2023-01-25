library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(ggrepel)
library(simplifyEnrichment)
library(ggalluvial)


#1 ORA
#@input:The list of gene symbols
#@output:The result of ORA for GO-BP
GO_enrich <- enrichGO(Genelist_GO, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod ='BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
GO_simplify <- simplify(GO_enrich,cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang")
GO_result <- as.data.frame(GO_simplify)



#2 GSEA
msigdbr_show_species()
Mm_terms = msigdbr(species = "Mus musculus")
Mm_arrange = Mm_terms %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)


GSEA_category=function(Genelist_GSEA,Categories=c("KEGG","REACTOME","CP","GO:BP","Hallmark")){
	#@Genelist_GSEA:The gene list of log2FoldChange whose name is gene symbols
	#@Categories:The categories selected for GSEA
	#@output:The result of GSEA
	if(Categories=="KEGG"){
		category_terms=msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
	}else if(Categories=="REACTOME"){
		category_terms=msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
	}else if(Categories=="CP"){
		category_terms=msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
	}else if(Categories=="GO:BP"){
		category_terms=msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
	}else if(Categories=="Hallmark"){
		category_terms=msigdbr(species = "Mus musculus", category = "H")
	}
	category_terms_select = category_terms %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
	Result_GSEA=GSEA(gene=Genelist_GSEA ,TERM2GENE = category_terms_select,pvalueCutoff = 0.05)
	return(Result_GSEA)
}

Bubbleplot_GSEA=function(Result_GSEA,Categories=c("KEGG","REACTOME","CP","GO:BP","Hallmark")){
	#@Result_GSEA:The result of GSEA
	#@Categories:The categories selected for GSEA
	#@output:The bubble plot of GSEA
	Result_GSEA=as.data.frame(Result_GSEA)
	Result=Result_GSEA[,c(1,5,8)]
	Result$NES=as.numeric(Result$NES)
	Result$qvalues=as.numeric(Result$qvalues)
	Result$log10=-log(Result$qvalues,10)
	Result$label=substr(Result$ID,6,nchar(Result$ID))
	Result$label=stringr::str_to_title(Result$label)
	Result$label=gsub("_"," ",Result$label)
	Bubble_P=ggplot(data=Result)+geom_point(aes(x=NES,y=log10,fill=NES,size=log10,alpha=log10),color="black",shape=21,show_guide=F)+scale_fill_gradientn(colours=c("#009ad6","#e0861a"))+
		geom_hline(yintercept=-log(0.05,10),color="#d71345",size=1.2,linetype="dashed")+
		geom_text_repel(data=Result[Result$log10 > -log(0.05,10),],aes(x=NES,y=log10,label=label),size=3)+
		theme_bw()+xlab("Nomalized Enrichment Score (NES)")+ylab("-log10 FDR q-val")+ggtitle(paste0("MsigDB ",Categories," gene sets"))+
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.text=element_text(size=rel(1.2)),axis.title=element_text(size=rel(1.2)),strip.text=element_text(size=rel(1.2)),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2)))
	return(Bubble_P)
}



#3 Sankey map of enrichment analysis

Enrich_cluster=function(GO_result_list,Order_TF){
	#@GO_result_list:The list R-object of results computed by clusterProfiler (like #1 ORA)
	#@Order_TF:The ordered TF(regulon) names (i.e. according to mad)
	#@output:The mapping of regulons, terms and clusters
	for(i in 1:length(GO_result_list)){
		name_i=names(GO_result_list)[i]
		GO_result_list[[i]]$TF_Name=rep(name_i,nrow(GO_result_list[[i]]))
	}
	GO_result_df=do.call(rbind,lapply(GO_result_list,function(x){return(x[,c(1,2,7,10)])}))
	GO_result_df=GO_result_df[,c(4,1,2)]
	GO_result_df=GO_result_df[GO_result_df$TF_Name %in% Order_TF,]
	Description=unique(GO_result_df$Description)
	GO_ID=unique(GO_result_df$ID)
	GO_simi=GO_similarity(GO_ID)
	GO_simp=simplifyGO(GO_simi)
	cluster=GO_simp$cluster
	for(j in 1:length(unique(cluster))){
		Des_j=Description[which(cluster==j)]
		Des_j=list(Des_j)
		names(Des_j)=paste0("cluster_",j)
		if(j==1){word_list=Des_j}else{word_list=c(word_list,Des_j)}
	}
	for(k in 1:length(word_list)){
		length_k=length(word_list[[k]])
		term_k=data.frame(cluster=rep(names(word_list)[k],length_k),terms=word_list[[k]])
		if(k==1){term_cluster=term_k}else{term_cluster=rbind(term_cluster,term_k)}
	}
	term_count=as.data.frame(table(term_cluster$cluster))
	term_lt_10=c()
	for(h in 1:length(word_list)){
		cluster_h=paste0("cluster_",h)
		number_h=term_count[term_count$Var1==cluster_h,2]
		if(number<=10){term_lt_10=c(term_lt_10,cluster_h)}
	}
	term_cluster$cluster_rename=term_cluster$cluster
	term_cluster$cluster_rename=ifelse(term_cluster$cluster %in% term_lt_10,"Others",term_cluster$cluster)
	GO_result_df$cluster=rep("Others",nrow(GO_result_df))
	GO_result_df$cluster=apply(GO_result_df,1,function(x){
		term_cluster[term_cluster$terms==x[3],3]
	})
	return(GO_result_df)
}


Sankey_map=function(GO_result_df,Order_TF,rename_cluster){
	#@GO_result_df:The output of 'Enrich_cluster' function
	#@Order_TF:The ordered TF(regulon) names (i.e. according to mad)
	#@rename_cluster:The new name of clusters
	#@output:The sankey map
	old_name=sort(unique(GO_result_df$cluster_rename))
	new_name=rename_cluster
	names(new_name)=old_name
	GO_result_df$cluster_rename <- new_name[as.character(old_name)]
	term_plot=GO_result_df[,c(1,4)]
	plot_df <- to_lodes_form(term_plot[,1:ncol(term_plot)],axes = 1:ncol(term_plot),id = "value")
	plot_df$stratum=factor(plot_df$stratum,levels=c(Order_TF,rename_cluster))
	p_sankey=ggplot(plot_df, aes(x = x, fill=stratum, label=stratum,stratum = stratum, alluvium  = value))+geom_flow(width = 0.3,curve_type = "sine",alpha = 0.5,color = 'white',size = 0.1)+geom_stratum(width = 0.28)+geom_text(stat = 'stratum', size = 2, color = 'black')+scale_fill_manual(values = col)+theme_void()+theme(legend.position = 'none')
	return(p_sankey)
}
