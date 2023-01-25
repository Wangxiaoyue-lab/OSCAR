library(Seurat)
library(scMAGeCK)

# 1 Prepare for scMAGeCK
##@input:The proceessed Seurat object

load(file="./pilot_DM_EM.Rdata")
DM2_filt=subset(DM_EM,subset=orig.ident=="DM")
save(DM2_filt,"./pilot_DM2_filt .RDS")
BARCODE <- DM2_filt@meta.data
BARCODE$cell <- rownames(BARCODE)
BARCODE_sub <- BARCODE[,c("cell","sgRNA_name","sgRNA_sequence","target_gene")]
colnames(BARCODE_sub) <- c("cell","barcode","sgrna","gene")
BARCODE_sub$read_count = rep(1,nrow(BARCODE_sub))
BARCODE_sub$umi_count = rep(1,nrow(BARCODE_sub))
BARCODE_sub$gene=BARCODE_sub$barcode
BARCODE_sub$gene[grep(BARCODE_sub$barcode,pattern="non_target")]="non_target"
BARCODE_sub$gene=gsub("-",".",BARCODE_sub$gene)
write.table(BARCODE_sub,file="./BARCODE_sub.txt",quote=F,sep="\t",col.names=T)



# 2 Run scMAGeCK
##@input:The proceessed Seurat object
##@output:The result of scMAGeCK

RDS="./pilot_DM2_filt.RDS"
BARCODE_3="./BARCODE_sub.txt"
RRAPATH <- NULL 
NEGCTRL="non_target"

lr_result <- scmageck_lr(BARCODE=BARCODE_3, RDS=RDS, LABEL='dox_scmageck_lr', NEGCTRL = NEGCTRL, PERMUTATION = 1000, SAVEPATH=NULL, LAMBDA=0.01,GENE_FRAC = 0)
lr_score <- lr_result[1][[1]] 
lr_score_pval <- lr_result[2][[1]]
save(lr_result,lr_score,lr_score_pval,file="./scmageck_lr.Rdata")

lr_score_sub=t(lr_score[,-1])
lr_score_sub=apply(lr_score_sub,1:2,as.numeric)
colnames(lr_score_sub)=stringr::str_to_title(colnames(lr_score_sub)) 
colnames(lr_score_sub)=sub("[\\.]","_",colnames(lr_score_sub))
colnames(lr_score_sub)=sub("NegCtrl","NTCs",colnames(lr_score_sub))
mad_order=apply(lr_score_sub,1,mad) #sd
score_order=lr_score_sub[order(mad_order,decreasing=T),]
score_order_top20=score_order[1:20,]
colnames(score_order_top20)[16]="NTCs"

