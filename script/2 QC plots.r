library(Seurat)
library(ggplot2)
library(patchwork)


# 1 The number of sgRNAs per cell of OSCAR
##@input:The counting table of sgRNAs per cell
sgrna_cell=read.csv("./cell_sgrna.csv",row.names=1,header=T)
sgrna_cell$sgrna_num=ifelse(sgrna_cell$sgrna_num>7,"≥8",as.character(sgrna_cell$sgrna_num))

#colour_fill=c(rep("#afb4db",6),rep("#76becc",3))
colour_fill=c(rep("#afb4db",1),rep("#76becc",1))

#for(i in 1:length(unique(sgrna_cell$sample))){
for(i in 1:length(unique(sgrna_cell$type))){
	#sample=unique(sgrna_cell$sample)[i]
	sample=unique(sgrna_cell$type)[i]
	#subset=sgrna_cell[sgrna_cell$sample==sample,]
	subset=sgrna_cell[sgrna_cell$type==sample,]
	subset$sgrna_num=factor(subset$sgrna_num,levels=c("0","1","2","3","4","5","6","7","≥8"))
	table_plot=as.data.frame(table(subset$sgrna_num))
	table_plot$rate=round(table_plot$Freq/sum(table_plot$Freq),4)
	p=ggplot(table_plot,aes(x=Var1,y=rate))+
		geom_bar(stat="identity",fill=colour_fill[i],color="black")+
		xlab("")+ylab(sample)+scale_y_continuous(breaks = c(0.0, 0.1, 0.2,0.3,0.4,0.5,0.6),limits = c(0.0,0.6))+
		theme_bw()+theme(plot.margin =unit(c(-0.75, 0, -0.75, 0), "cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.text.x=element_blank(),axis.text.y=element_text(size=5),axis.ticks.x=element_blank(),axis.title.y=element_text(size=10, angle=0, vjust=0.5))
 
	plot_l=list(p)
	if(i==1){plot_list=plot_l}else{plot_list=c(plot_list,plot_l)}
}
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +theme(axis.text.x=element_text(size=10, vjust = 0.5, hjust = 0.5), axis.ticks.x = element_line())
 
pdf(file="./DM+EM_type_sgrna_per_cell.pdf",width=6,height=4)
q <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
q
dev.off()





# 2 The number of cells per sgrna
##@input:The proceesed Seurat object
load(file="./filt_dm_em2.Rdata")

unique_meta=data.filt@meta.data
unique_meta_DM=unique_meta[unique_meta$type=="DM",]

DM_gene=as.data.frame(table(unique_meta_DM$gene))
DM_gene$log2Freq=log(DM_gene$Freq,2)

pdf("./DM_gene_histogram.pdf",width=4,height=4)
ggplot(DM_gene,aes(x=log2Freq))+geom_histogram(bins=60,fill="grey",colour="black")+geom_density(size=1,colour="green")+
	xlab("Log2")+ylab("Count")+scale_x_continuous(limits=c(0,20))+ggtitle("Cell number per gene target")+
	theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white",color="black"))
dev.off()


DM_sgrna=as.data.frame(table(unique_meta_DM$barcode))
DM_sgrna$log2Freq=log(DM_sgrna$Freq,2)

pdf("./DM_sgrna_histogram.pdf",width=4,height=4)
ggplot(DM_sgrna,aes(x=log2Freq))+geom_histogram(bins=60,fill="grey",colour="black")+geom_density(size=1,colour="green")+
	xlab("Log2")+ylab("Count")+scale_x_continuous(limits=c(0,20))+ggtitle("Cell number per sgRNA target")+
	theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white",color="black"))
dev.off()



# 3 The bubble plot
##@input:The proceesed Seurat object
load(file="./filt_dm_em2.Rdata")

unique_meta=data.filt@meta.data[,c(2,4,10,11)]
bubble=as.data.frame(table(unique_meta$type,unique_meta$gene))
colnames(bubble)=c("type","gene","num")
bubble$type_value=ifelse(bubble$type=="DM",0.35,0.6)
bubble$type_value=as.numeric(bubble$type_value)

bubble_DM=bubble[bubble$type=="DM",]
bubble_DM$order=1:nrow(bubble_DM)
bubble_EM=bubble[bubble$type=="EM",]
bubble_EM$order=1:nrow(bubble_EM)
Angle_gt_90=round(nrow(bubble_DM)*90/360,0)*(360/nrow(bubble_DM))
bubble_DM$angle=apply(bubble_DM,1,function(x){Angle_gt_90-as.numeric(x[5])*360/nrow(bubble_DM) })
bubble_EM$angle=apply(bubble_EM,1,function(x){Angle_gt_90-as.numeric(x[5])*360/nrow(bubble_DM)})
bubble_DM$angle=ifelse(bubble_DM$angle< -90,bubble_DM$angle+540,bubble_DM$angle)
bubble_EM$angle=ifelse(bubble_EM$angle< -90, bubble_EM$angle+540,bubble_EM$angle)
bubble=rbind(bubble_DM,bubble_EM)
bubble$hjust=ifelse(bubble$angle > 270 &  bubble$angle < 450, 1,0)
 
pdf("./bubble.pdf",width=10,height=10)
ggplot(data=bubble)+geom_segment(aes(x=gene,xend=gene),y=-0.6,yend=0.75,size=0.5,colour="#d3d7d4",alpha=0.4)+
	geom_point(aes(size=num,x=gene,y=type_value,colour=factor(type)))+scale_colour_manual(values=c("DM"="#50b7c1","EM"="#afb4db"))+
	geom_text(aes(x=gene,label=num,y=type_value+0.03,angle=angle,hjust=hjust),color="#130c0e")+
	geom_text(aes(x=gene,label=gene,angle=angle,hjust=hjust),color="#130c0e",y=0.03)+
	geom_hline(yintercept=c(-0.2,0,0.3,0.55,0.75),color="#d3d7d4",size=1)+
	xlab("")+ylab("")+ylim(-0.6,0.75)+
	theme_bw()+theme(panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour="black",size=0),axis.ticks =element_blank(),axis.text = element_blank())+coord_polar(theta = "x")
dev.off()
