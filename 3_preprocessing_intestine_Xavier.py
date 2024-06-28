### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_prefix> <out_dir>\n"
arg1 = "\nin_dir: input dir containing the genes\n"
arg2 = "\nin_obj: input h5ad file\n"

parser = OptionParser(usage=usage + arg1 + arg2)
(options, args) = parser.parse_args()

if len(args) < 2:
	sys.stderr.write("ERROR: <in_prefix> and <out_dir> are required\nTry --help for help\n")
	exit()
in_prefix = args[0]
out_dir = args[1]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * # TODO: import from the (external) git directory

for cond in ['CO', 'TI']:
    in_obj = in_prefix + "_" + cond + ".h5ad"
    adata = preprocess_adata_intestine_Xavier(in_obj, cond = cond)
    out_obj = out_dir + "/Intestine_Xavier_" + cond + ".h5ad"
    adata.write_h5ad(out_obj)

exit()

