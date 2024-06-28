### conda activate scvi-env

### options

EXP = "Smillie et al. 2019"

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <in_obj> <metadata> <cells>\n"
arg1 = "\nin_obj: input h5ad file containing the atlas\n"
arg2 = "\nmetadata: input data frame containing the metadata, including cell name and donor_id\n"
arg3 = "\ncells: output containing the list of cells to retain\n"

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_obj>, <metadata> and <cells> are required\nTry --help for help\n")
	exit()
in_obj = args[0]
in_dir = args[1]
out = args[2]

### environment

import numpy as np
import scanpy as sc
import pandas as pd
from src.cross_tissue_fibrosis.datasets import * 

donors = extract_donor_from_adata_cellTypist(in_obj, exp = EXP)
cells = extract_cells_intestine_Xavier(in_dir, to_exclude = donors)
df = pd.DataFrame(cells)
df.to_csv(out, sep = "\t")

exit()

