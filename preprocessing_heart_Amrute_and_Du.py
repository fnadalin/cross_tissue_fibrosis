### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_dir_1> <in_dir_2> <out_obj>\n"
arg1 = "\nin_dir_1: directory containing input files (dataset 1)\n"
arg2 = "\nin_dir_2: directory containing input files (dataset 2)\n"
arg3 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_dir_1>, <in_dir_2> and <out_obj> are required\nTry --help for help\n")
	exit()
in_dir_1 = args[0]
in_dir_2 = args[1]
out_obj = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * 

adata = preprocess_adata_heart_Amrute_filt_by_meta(in_dir_1)
adata_2 = preprocess_adata_heart_Du(in_dir_2)
adata = ad.concat([adata, adata_2], merge = "same")
adata.write_h5ad(out_obj)

exit()

