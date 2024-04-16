### conda activate scvi-env

"""
Learn the latent space on the reference (unsupervised or supervised)
Map the query onto the reference
Predict cell types (only in the supervised version)
TODO: deal with soft labels in the reference (for cross-tissue cell-type prediction in the reference)
"""

def param_val(param_file, val):
    v = params[0] == val
    return params.loc[v][1][np.where(v)[0][0]]

### options

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

print("==== Latent embedding ====")

### environment

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from src.cross_tissue_fibrosis.anndata_objects import *
from src.cross_tissue_fibrosis.latent_space_learning import *
from src.cross_tissue_fibrosis.paths import *
# from src.cross_tissue_fibrosis.differential_abundance import *

### parameters
params = pd.read_csv(PARAMS, sep = "=", header = None)
folder      = param_val(params, "folder")
r_obj       = param_val(params, "r_obj")
q_obj       = param_val(params, "q_obj")
r_name      = param_val(params, "r_name")
q_name      = param_val(params, "q_name")
r_batch_key = param_val(params, "r_batch_key")
q_batch_key = param_val(params, "q_batch_key")
n_top_genes = int(param_val(params, "n_top_genes"))
flavor      = param_val(params, "flavor")
meta        = param_val(params, "meta").split(",")

# print params
params.to_csv(os.path.join(folder, "params.tsv"))

# paths to model dirs
r_mod_dir = create_mod_path(folder, r_name, r_batch_key, n_top_genes, flavor)
q_mod_dir = create_mod_path(folder, q_name, q_batch_key, n_top_genes, flavor)

# paths to object dirs
r_obj_dir = create_obj_path(folder, r_name, r_batch_key, n_top_genes, flavor)
q_obj_dir = create_obj_path(folder, q_name, q_batch_key, n_top_genes, flavor)

# preprocessing
print("Prepare anndata objects")
adata_ref = prepare_adata_for_scVI(r_obj, r_batch_key, n_top_genes = n_top_genes, flavor = flavor)
adata_query = prepare_adata_for_scVI(q_obj, r_batch_key, compute_hvg = False)
# hvg_r = ENSEMBL_features_from_adata(adata_ref) # CellTypist objects do not contain ENSEMBL IDs!!!
hvg_r = features_from_adata(adata_ref)
# adata_query = subset_by_ENSEMBL_features(adata_query, hvg_r)
adata_query = subset_by_features(adata_query, hvg_r) # now the query contains only a subset of hvg_r
# Add the missing features to the query, otherwise scvi gives an error!!!
# note: this is done after subsetting by hvg_r, so I don't have to create a huge temporary anndata object 
adata_merged = ad.concat([adata_ref, adata_query], join="outer")
obs_names = adata_query.obs.columns
# uns_names = adata_query.uns
adata_query = adata_merged[adata_query.obs.index,:].copy() # now the query contains all hvg_r, with padding zeroes for those hvgs only in ref (N.B. if the reference does not have ENSEMBL gene IDs, this information will be entirely lost!!)
del adata_merged
adata_query.obs = adata_query.obs[obs_names]
# adata_query.uns = adata_query.uns[uns_names] # FIXME: this throws an error

###### STRATEGY 1: LEARN LATENT SPACE ON QUERY #######

# learn scVI latent space on query (use the query as a reference)
# FIXME: this implies query-specific feature selection, so it must be saved as a separate object
# q_ref_mod = os.path.join(q_mod_dir, "scvi")
# run_scVI_reference(adata_query, q_ref_mod)

###### STRATEGY 2: LEARN LATENT SPACE ON REFERENCE #######

# learn scVI latent space on reference
print("learn scVI latent space on reference [scVI]")
r_mod = create_scvi_mod_ref(r_mod_dir)
run_scVI_reference(adata_ref, r_mod)

# update with query
print("update with query [scVI]")
q_mod = create_scvi_mod_query(q_mod_dir, r_name)
run_scVI_query(adata_query, r_mod, q_mod, r_name)

###### STRATEGY 3: LEARN LATENT SPACE ON REFERENCE + COVARIATE #######

for m in meta:
    # learn scANVI latent space on reference
    print("learn scVI latent space on reference [scANVI " + m + "]")
    r_an_mod = create_scanvi_mod_ref(r_mod_dir, m)
    run_scANVI_reference(adata_ref, r_mod, r_an_mod, m)

    # update with query 
    print("update with query [scANVI " + m + "]")
    q_an_mod = create_scanvi_mod_query(q_mod_dir, r_name, m)
    run_scANVI_query(adata_query, r_an_mod, q_an_mod, m, r_name)

    # predict cell type labels
    print("predict cell type labels [scANVI " + m + "]")
    labels_hard, labels_soft = predict_scANVI_query_labels(adata_query, q_an_mod)

    # paths to label files
    q_lab_hard, q_lab_soft = create_labels_file(q_obj_dir, r_name, m)

    # write labels to file
    pd.DataFrame(labels_hard).to_csv(q_lab_hard, sep = "\t")
    labels_soft.to_csv(q_lab_soft, sep = "\t")
    
###### SAVE LATENT EMBEDDINGS #######

# update objects with scVI latent embeddings
print("save embeddings")
# adata_query = update_adata_scVI_reference(adata_query, q_ref_mod)
adata_ref = update_adata_scVI_reference(adata_ref, r_mod)
adata_query = update_adata_scVI_query(adata_query, q_mod, r_name)

# update objects with scANVI latent embeddings
for m in meta:
    adata_ref = update_adata_scANVI_reference(adata_ref, r_an_mod, m)
    adata_query = update_adata_scANVI_query(adata_query, q_an_mod, r_name, m)

# paths to objects
r_obj = create_obj_file(r_obj_dir)
q_obj = create_obj_file(q_obj_dir)

# write objects to file
adata_ref.write_h5ad(r_obj)
adata_query.write_h5ad(q_obj)

exit()


