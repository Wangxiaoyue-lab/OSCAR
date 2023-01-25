##@tfs   mm_mgi_tfs.txt
##@feather  mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
##@tbl  motifs-v9nr_clust-nr.hgnc-m0.001-o0.0.tbl 
##@input_loom  data.filt.sample.loom

ls $tfs  $feather  $tbl  

pyscenic grn \
--num_workers 20 \
--output adj.data.filt.sample.tsv \
--method grnboost2 \
$input_loom  $tfs 


pyscenic ctx \
adj.data.filt.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output data.filt.reg.csv \
--num_workers 20  \
--mask_dropouts


pyscenic aucell \
$input_loom \
data.filt.reg.csv \
--output data.filt.out_SCENIC.loom \
--num_workers 20 