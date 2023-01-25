# transform csv to loom
##@input:The matrix obtained from Seurat object
##@output:The loom file to pySCENIC
import os, sys
os.getcwd() 
os.listdir(os.getcwd())  

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("data.filt.sample.csv"); 
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("data.filt.sample.loom",x.X.transpose(),row_attrs,col_attrs);