### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_obj> <category_file> <out_obj>\n"
arg1 = "\nin_obj: input h5ad file containing the original anndata file from CellTypist db\n"
arg2 = "\ncategory_file: tsv file containing the mapping between cell types and broad category\n"
arg3 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_obj>, <category_file> and <out_obj> are required\nTry --help for help\n")
	exit()
in_obj = args[0]
category_file = args[1]
out_obj = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * 

adata = preprocess_adata_cellTypist(in_obj, category_file)
adata.write_h5ad(out_obj)

exit()

