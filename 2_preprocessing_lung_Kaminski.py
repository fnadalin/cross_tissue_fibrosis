### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_dir> <in_obj> <out_obj>\n"
arg1 = "\nin_dir: input dir containing the genes\n"
arg2 = "\nin_obj: input h5ad file\n"
arg3 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_obj> and <out_obj> are required\nTry --help for help\n")
	exit()
in_dir = args[0]
in_obj = args[1]
out_obj = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * # TODO: import from the (external) git directory

adata = preprocess_adata_lung_Kaminski(in_dir, in_obj)
adata.write_h5ad(out_obj)

exit()

