### conda activate sklearn-env

"""
Predict DA and DE between normal and disease as a linear combination of cell type annotations in neighbourhoods
"""

def param_val(param_file, val):
    v = params[0] == val
    return params.loc[v][1][np.where(v)[0][0]]

### options

C_vec = [1,0.5,0.1,0.05,0.01,0.005,0.001] # inverse of the regularization weight

from optparse import OptionParser
import os
import sys

usage = "python %prog [options] <params>\n"
arg1 = "\nparams: file containing parameter values\n"

parser = OptionParser(usage=usage + arg1)
(options, args) = parser.parse_args()

if len(args) < 1:
	sys.stderr.write("ERROR: <params> is required\nTry --help for help\n")
	exit()
PARAMS=args[0]

### environment

import os
import numpy as np
import pandas as pd
from src.cross_tissue_fibrosis.multinomial_regression import *
from src.cross_tissue_fibrosis.paths import *

### parameters
params = pd.read_csv(PARAMS, sep = "=", header = None)
folder      = param_val(params, "folder")
r_name      = param_val(params, "r_name")
q_name      = param_val(params, "q_name")
q_batch_key = param_val(params, "q_batch_key")
n_top_genes = int(param_val(params, "n_top_genes"))
flavor      = param_val(params, "flavor")
meta        = param_val(params, "meta").split(",")

# paths to object dirs
q_milo_dir = create_milo_path(folder, q_name, q_batch_key, n_top_genes, flavor)

###### STRATEGY 1: LEARN LATENT SPACE ON QUERY #######

# TODO

###### STRATEGY 2: LEARN LATENT SPACE ON REFERENCE #######

q_milo = create_scvi_mod_query(q_milo_dir, r_name)
out_miloR_scvi_prefix = q_milo + "_miloR"
out_miloDE_scvi_prefix = q_milo + "_miloDE"
for m in meta:
    da_results_ann_file = out_miloR_scvi_prefix + "_" + r_name + "_" + m + ".tsv"
    de_results_ann_file = out_miloDE_scvi_prefix + "_" + r_name + "_" + m + ".tsv"
    for C in C_vec:
        coef, names = da_logit_from_celltype(da_results_ann_file, C)
        write_da_logit_from_celltype(out_miloR_scvi_prefix + "_" + r_name + "_" + m, coef, names, C)
        coef_df = de_logit_from_celltype(de_results_ann_file, C)
        write_de_logit_from_celltype(out_miloDE_scvi_prefix + "_" + r_name + "_" + m, coef_df, C)

###### STRATEGY 3: LEARN LATENT SPACE ON REFERENCE + COVARIATE #######

for m in meta:
    q_milo = create_scanvi_mod_query(q_milo_dir, r_name, m) 
    out_miloR_scanvi_prefix = q_milo + "_miloR"
    out_miloDE_scanvi_prefix = q_milo + "_miloDE"
    da_results_ann_file = out_miloR_scanvi_prefix + "_" + r_name + "_" + m + ".tsv"
    de_results_ann_file = out_miloDE_scanvi_prefix + "_" + r_name + "_" + m + ".tsv"
    for C in C_vec:
        coef, names = da_logit_from_celltype(da_results_ann_file, C)
        write_da_logit_from_celltype(out_miloR_scanvi_prefix + "_" + r_name + "_" + m, coef, names, C)
        coef_df = de_logit_from_celltype(de_results_ann_file, C)
        write_de_logit_from_celltype(out_miloDE_scanvi_prefix + "_" + r_name + "_" + m, coef_df, C)

exit()




