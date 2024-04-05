### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <cell_obj> <nuclei_obj> <out_obj>\n"
arg1 = "\ncell_obj: input h5ad file containing cells\n"
arg2 = "\nnuclei_obj: input h5ad file containing nuclei\n"
arg3 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <cell_obj>, <nuclei_obj> and <out_obj> are required\nTry --help for help\n")
	exit()
cell_obj = args[0]
nuclei_obj = args[1]
out_obj = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * # TODO: import from the (external) git directory

adata = preprocess_adata_heart_Lavine(cell_obj, nuclei_obj)
adata.write_h5ad(out_obj)

exit()

