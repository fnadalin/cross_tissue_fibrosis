# N.B.: obs rownames are removed upon conversion to R adata!!! Save the barcodes as an adata.obs field.
# N.B.: better to set the feature names as ENSEMBL instead of SYMBOL (even though they are both saved in adata.var)

### environment

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import warnings
import re
import os
# from utils import pandas_rename_columns

### global settings

metadata = ['dataset', 'sample_id', 'cell_type', 'condition', 'tissue', 'category']
celltypist_fields = ['Dataset', 'donor_id', 'Curated_annotation', 'disease', 'tissue', 'category']
Lavine_fields = ['dataset', 'orig.ident', 'Names', 'Condition', 'tissue', 'category']
Kaminski_fields = ['dataset', 'Subject_Identity', 'Subclass_Cell_Identity', 'Disease_Identity', 'tissue', 'CellType_Category']

Henderson_healthy = ["GSM4041150_healthy1_cd45+", 
                     "GSM4041151_healthy1_cd45-A", 
                     "GSM4041152_healthy1_cd45-B", 
                     "GSM4041153_healthy2_cd45+", 
                     "GSM4041154_healthy2_cd45-", 
                     "GSM4041155_healthy3_cd45+", 
                     "GSM4041156_healthy3_cd45-A", 
                     "GSM4041157_healthy3_cd45-B", 
                     "GSM4041158_healthy4_cd45+", 
                     "GSM4041159_healthy4_cd45-", 
                     "GSM4041160_healthy5_cd45+"] 
Henderson_cirrhotic = ["GSM4041161_cirrhotic1_cd45+", 
                       "GSM4041162_cirrhotic1_cd45-A", 
                       "GSM4041163_cirrhotic1_cd45-B", 
                       "GSM4041164_cirrhotic2_cd45+", 
                       "GSM4041165_cirrhotic2_cd45-", 
                       "GSM4041166_cirrhotic3_cd45+", 
                       "GSM4041167_cirrhotic3_cd45-", 
                       "GSM4041168_cirrhotic4_cd45+", 
                       "GSM4041169_cirrhotic5_cd45+"]
                       
Kaminski_genes = "GSE136831_AllCells.GeneIDs.txt"

### functions

def preprocess_adata_cellTypist(in_obj, category_file):
    """
    Put the cellTypist adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    obs - retain only "sample_ID", "condition", "dataset"
    Cell ID must be unique 
    
    Parameters
    ----------
    in_obj
        adata object filename
    category_file
        tsv file containing the cell-type on column 1 and the category on column 2
    """
    adata = sc.read_h5ad(in_obj)
    map_ct_cat = pd.read_csv(category_file, sep = "\t", header = None, index_col = False)
    
    # extract and rename obs fields
    idx = R_match(adata.obs.columns, np.array(celltypist_fields))
    df = adata.obs.loc[:,idx]
    adata.obs = pandas_rename_columns(df, celltypist_fields, metadata)
    
    # check that all cell types in the object are represented in the cell-type/category table
    ct_unique = pd.unique(adata.obs.loc[:,'cell_type'])
    matched = sum(map_ct_cat.iloc[i,0] in ct_unique for i in range(map_ct_cat.shape[0]))
    if matched < len(ct_unique):
        warnings.warn(str(len(ct_unique)-matched) + " cell-types will not be assigned a category")
    
    # map cell-types to categories
    v = adata.obs.loc[:,'cell_type']
    w = map_ct_cat.iloc[:,0]
    idx = R_match_idx(v.to_numpy(), w.to_numpy())
    adata.obs.loc[:,'category'] = R_match(map_ct_cat.iloc[idx,1], v.index)
    
    return adata


# GSE183852
def preprocess_adata_heart_Lavine(cell_obj, nuclei_obj):
    """
    Put the Lavine adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    obs - retain only "sample_ID", "condition", "dataset"
    Cell ID must be unique 
    
    Parameters
    ----------
    cell_obj
        adata object filename containing cells
    nuclei_obj
        adata object filename containing nuclei
    """
    adata = sc.read_h5ad(cell_obj)
    adata_nuclei = sc.read_h5ad(nuclei_obj)
    adata = ad.concat([adata, adata_nuclei])
    
    # extract, integrate and rename obs fields
    df = adata.obs.loc[:,Lavine_fields]
    adata.obs['dataset'] = 'Lavine'
    adata.obs['tissue'] = 'heart'
    adata.obs['category'] = 'Unknown'
    adata.obs = pandas_rename_columns(df, Lavine_fields, metadata)
    
    return adata.copy()


# GSE136103
# N.B. Input files have been gunzipped so that python recognises them as "legacy" cellRanger output
def preprocess_adata_liver_Henderson(in_dir):
    """
    Put the Henderson adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    obs - retain only "sample_ID", "condition", "dataset"
    Cell ID must be unique 
    
    Parameters
    ----------
    in_dir
        input folder containing the matrices (*genes*, *barcodes*, *matrix*)
    """
    # create adata
    adata = concat_adata_from_multi_mtx(in_dir, Henderson_healthy)
    adata.obs['condition'] = 'normal'
    adata_cirrhotic = concat_adata_from_multi_mtx(in_dir, Henderson_cirrhotic)
    adata_cirrhotic.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_cirrhotic], merge = "same")
    adata.obs['dataset'] = 'Henderson'
    adata.obs['tissue'] = 'liver'
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['category'] = 'Unknown'
    
    return adata.copy()


# GSE136831
# matrix preprocessed with R function load_Kaminski() because it's faster than python ad.read_mtx
def preprocess_adata_lung_Kaminski(in_dir, in_obj):
    """
    Put the Kaminski adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    obs - retain only "sample_ID", "condition", "dataset"
    Cell ID must be unique 
    
    Parameters
    ----------
    in_dir
        input directory containing genes
    in_obj
        input adata object obtained from Seurat
    """
    genes_file = os.path.join(in_dir, Kaminski_genes)
    genes = pd.read_csv(genes_file, sep = "\t", header = 0, index_col = False)
    genes.index = genes['HGNC_EnsemblAlt_GeneID']
    genes.index.name = None

    # process features
    adata = sc.read_h5ad(in_obj)
    adata.var.index = adata.var['features']
    adata.var.index.name = None
    adata.var["gene_ids"] = genes["Ensembl_GeneID"]

    # process cells and metadata
    adata.obs['dataset'] = 'Kaminski'
    adata.obs['tissue'] = 'lung'
    adata.obs = pandas_rename_columns(adata.obs, Kaminski_fields, metadata)

    # extract cells belonging to selected conditions
    adata = adata[adata.obs['condition'].isin(['Control', 'IPF'])].copy()
    idx_normal = np.where(adata.obs['condition'] == 'Control')
    idx_disease = np.where(adata.obs['condition'] == 'IPF')
    adata.obs['condition'].iloc[idx_normal] = 'normal'
    adata.obs['condition'].iloc[idx_disease] = 'disease'
    
    return adata.copy()
    

def concat_adata_from_multi_mtx(in_dir, prefixes):
    """
    Create an anndata object from multiple mtx files
    
    Parameters
    ----------
    in_dir
        input folder containing the matrices (*genes*, *barcodes*, *matrix*)
    prefixes
        list of file prefixes, one per matrix
    """
    init = False
    for prefix in prefixes:
        adata_tmp = sc.read_10x_mtx(in_dir, make_unique = True, prefix = prefix + '_')
        s = re.sub(r'GSM[0-9]+_', '', prefix)
        sample_name = re.sub(r'_.+', '', s)
        exp_name = re.sub(r'_.+', '', prefix)
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index
        adata_tmp.obs['library'] = prefix
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()



