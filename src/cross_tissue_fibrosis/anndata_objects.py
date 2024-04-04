### environment

import scanpy as sc

### global settings

# batch_key = "Diagnosis"

### functions

def prepare_adata_for_scVI(in_obj, batch_key, compute_hvg = True, n_top_genes = 5000, flavor = "seurat_v3"):    
    """
    Compute hvgs and return the subsetted object
    
    Parameters
    ----------
    in_obj
        adata object filename
    batch_key
        compute hvgs within each batch_key and then take the union across batches
    """    
    adata = sc.read_h5ad(in_obj)
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.raw = adata
    if compute_hvg:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes = n_top_genes, 
            subset = True,
            layer = "counts",
            flavor = flavor,
            batch_key = batch_key
        )
    return adata.copy()


def subset_by_obs(adata, meta, values):    
    """
    Subset the object to the cells that match values in meta
    
    Parameters
    ----------
    adata
        input adata object 
    meta
        field in adata.obs
    values
        list of values in "meta" to subset on
    """        
    index = adata.obs[meta].isin(values)
    cells = adata.obs.index[index]
    adata = adata[cells,:]
    return adata.copy()


def subset_by_features(adata, features):
    adata = adata[:,features]
    return adata.copy()


def features_from_adata(adata):
    return adata.var.index


def subset_by_ENSEMBL_features(adata, features):
    index = adata.var["gene_ids"].isin(features)
    features = adata.var.index[index]
    adata = adata[:,features]
    return adata.copy()


def ENSEMBL_features_from_adata(adata):
    return adata.var["gene_ids"]


