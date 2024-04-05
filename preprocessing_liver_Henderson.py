### conda activate scvi-env

### options

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_dir> <out_obj>\n"
arg1 = "\nin_dir: directory containing input files\n"
arg2 = "\nout_obj: output h5ad file formatted for our pipeline\n"

parser = OptionParser(usage=usage + arg1 + arg2)
(options, args) = parser.parse_args()

if len(args) < 2:
	sys.stderr.write("ERROR: <in_dir> and <out_obj> are required\nTry --help for help\n")
	exit()
in_dir = args[0]
out_obj = args[1]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * 

adata = preprocess_adata_liver_Henderson(in_dir)
adata.write_h5ad(out_obj)

exit()

