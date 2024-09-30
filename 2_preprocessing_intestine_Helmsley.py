### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <normal_obj> <disease_obj> <out_obj>\n"
arg1 = "\nnormal_obj: input h5ad file containing normal samples\n"
arg2 = "\ndisease_obj: input h5ad file containing disease samples\n"
arg3 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <normal_obj>, <disease_obj> and <out_obj> are required\nTry --help for help\n")
	exit()
normal_obj = args[0]
disease_obj = args[1]
out_obj = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * # TODO: import from the (external) git directory

adata = preprocess_adata_intestine_Helmsley(normal_obj, disease_obj)
adata.write_h5ad(out_obj)

exit()

