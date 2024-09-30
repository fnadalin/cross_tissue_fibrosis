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
from src.cross_tissue_fibrosis.utils import *

### global settings

metadata = ['dataset', 'sample_id', 'cell_type', 'condition', 'tissue', 'category']
celltypist_fields = ['Dataset', 'donor_id', 'Curated_annotation', 'disease', 'tissue', 'category']
Lavine_fields = ['dataset', 'orig.ident', 'Names', 'Condition', 'tissue', 'category']
Kaminski_fields = ['dataset', 'Subject_Identity', 'Subclass_Cell_Identity', 'Disease_Identity', 'tissue', 'CellType_Category']
Xavier_fields = ['dataset', 'biosample_id', 'Celltype', 'disease', 'tissue', 'category']
Helmsley_fields = ['dataset', 'sample.name',  'cell_type_final', 'condition', 'tissue', 'category']

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
                       
Sethupathy_healthy = [ "GSM5024090_nIBD_17_0236_206",
                       "GSM5024091_nIBD_17_0236_214",
                       "GSM5024092_nIBD_17_0236_216",
                       "GSM5024093_nIBD_17_0236_217" ]
Sethupathy_cd = [ "GSM5024087_CD_uninflamed_17_0236_189",
                  "GSM5024088_CD_uninflamed_17_0236_299",
                  "GSM5024089_CD_uninflamed_17_0236_364" ]
                       
Kaminski_genes = "GSE136831_AllCells.GeneIDs.txt"

Xavier_meta = "scp_metadata_combined.v2.txt"
Xavier_genes_CO = "CO_EPI.scp.features.tsv"
Xavier_genes_TI = "TI_EPI.scp.features.tsv"

# removed corrupted file GSM8352066_MA29
Amrute_healthy = [ "GSM8352048_MA5",
                   "GSM8352049_MA6",
                   "GSM8352055_MA13",
                   "GSM8352060_MA23",
                   "GSM8352063_MA26",
                   "GSM8352070_MA33" ]

Amrute_disease = [ "GSM8352050_MA7",
                   "GSM8352051_MA8",
                   "GSM8352052_MA9",
                   "GSM8352053_MA10",
                   "GSM8352054_MA11",
                   "GSM8352056_MA14",
                   "GSM8352057_MA19",
                   "GSM8352058_MA20",
                   "GSM8352059_MA22",
                   "GSM8352061_MA24",
                   "GSM8352062_MA25",
                   "GSM8352064_MA27",
                   "GSM8352065_MA28",
                   "GSM8352067_MA30",
                   "GSM8352068_MA31",
                   "GSM8352069_MA32" ]

Amrute_metadata = "GSE270788_metadata.csv"

Du_healthy = [ "Control_1_filtered_feature_bc_matrix",
               "Control_2_filtered_feature_bc_matrix" ]

Du_disease = [ "Heart_failure_1_filtered_feature_bc_matrix",
               "Heart_failure_2_filtered_feature_bc_matrix",
               "Heart_failure_3_filtered_feature_bc_matrix",
               "Heart_failure_4_filtered_feature_bc_matrix" ]

### functions

def preprocess_adata_cellTypist(in_obj, category_file, datasets = None):
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


def extract_donor_from_adata_cellTypist(in_obj, exp):
    adata = sc.read_h5ad(in_obj)
    v = (adata.obs['Dataset'] == exp)
    donors = np.unique(adata.obs.iloc[np.where(v)]['donor_id'])
    
    return donors


def extract_cells_intestine_Xavier(in_dir, to_exclude = [], to_keep = []):
    if to_keep == "" and to_exclude == "":
        raise ValueError("Either to_keep or to_exclude must be non-empty strings")
    df = in_dir + "/" + Xavier_meta
    data = pd.read_csv(df, sep = "\t", header = 0, dtype = "str")
    data = data.iloc[1:data.shape[0],:]
#    data['donor_id'] = data['donor_id'].astype("str")
    
    to_exclude = np.asarray(to_exclude).astype("str")
    to_keep = np.asarray(to_keep).astype("str")
    
    v = np.ones(data.shape[0], dtype = "bool")
    if len(to_keep) > 0:
        v = v == (data['donor_id'].isin(to_keep).to_numpy())
    if len(to_exclude) > 0:
        v = np.not_equal(v, data['donor_id'].isin(to_exclude).to_numpy())
    cells = data.iloc[np.where(v)]['NAME'].to_numpy() 
    
    return cells


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
    adata.obs['dataset'] = 'Lavine'
    adata.obs['tissue'] = 'heart'
    adata.obs['category'] = 'Unknown'
    df = adata.obs.loc[:,Lavine_fields]
    adata.obs = pandas_rename_columns(df, Lavine_fields, metadata)
    adata.obs['cell_type'] = adata.obs['cell_type'].to_numpy().astype('str')
    
    # rename conditions
    di = {'Donor': 'normal', 'DCM': 'disease'}
    adata.obs["condition"].replace(di, inplace=True)
    
    del adata.raw # without this, LoadH5Seurat from the h5Seurat object obtained from h5ad will fail!!!
    return adata.copy()


def preprocess_adata_intestine_Helmsley(normal_obj, disease_obj):
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
    adata_n = sc.read_h5ad(normal_obj)
    adata_d = sc.read_h5ad(disease_obj)
    adata = ad.concat([adata_n, adata_d])
    
    # extract, integrate and rename obs fields
    adata.obs['dataset'] = 'Helmsley'
    adata.obs['tissue'] = 'intestine'
    adata.obs['category'] = 'Unknown'
    df = adata.obs.loc[:,Helmsley_fields]
    adata.obs = pandas_rename_columns(df, Helmsley_fields, metadata)
    adata.obs['cell_type'] = adata.obs['cell_type'].to_numpy().astype('str')

    # rename conditions
    di = {'control': 'normal', 'disease': 'disease'}
    adata.obs["condition"].replace(di, inplace=True)
    
    del adata.raw # without this, LoadH5Seurat from the h5Seurat object obtained from h5ad will fail!!!
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
    adata = concat_adata_from_multi_mtx_liver_Henderson(in_dir, Henderson_healthy)
    adata.obs['condition'] = 'normal'
    adata_cirrhotic = concat_adata_from_multi_mtx_liver_Henderson(in_dir, Henderson_cirrhotic)
    adata_cirrhotic.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_cirrhotic], merge = "same")
    adata.obs['dataset'] = 'Henderson'
    adata.obs['tissue'] = 'liver'
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['category'] = 'Unknown'
    
    return adata.copy()


# GSE247468
def preprocess_adata_heart_Du(in_dir):
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
    adata = concat_adata_from_multi_mtx_heart_Du(in_dir, Du_healthy)
    adata.obs['condition'] = 'normal'
    adata_disease = concat_adata_from_multi_mtx_heart_Du(in_dir, Du_disease)
    adata_disease.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_disease], merge = "same")
    adata.obs['dataset'] = 'Du'
    adata.obs['tissue'] = 'heart'
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['category'] = 'Unknown'
    
    return adata.copy()
    

# GSE164985
def preprocess_adata_intestine_Sethupathy(in_dir):
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
    adata = concat_adata_from_multi_mtx_intestine_Sethupathy(in_dir, Sethupathy_healthy)
    adata.obs['condition'] = 'normal'
    adata_cd = concat_adata_from_multi_mtx_intestine_Sethupathy(in_dir, Sethupathy_cd)
    adata_cd.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_cd], merge = "same")
    adata.obs['dataset'] = 'Sethupathy'
    adata.obs['tissue'] = 'intestine'
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
    df = adata.obs.loc[:,Kaminski_fields]
    adata.obs = pandas_rename_columns(df, Kaminski_fields, metadata)
    
    # extract cells belonging to selected conditions
    adata = adata[adata.obs['condition'].isin(['Control', 'IPF'])].copy()
    
    # rename conditions
    di = {'Control': 'normal', 'IPF': 'disease'}
    adata.obs["condition"].replace(di, inplace=True)
    
    del adata.raw # this solves the bug in adata.write_h5ad(): "ValueError: '_index' is a reserved name for dataframe columns."
    return adata.copy()


# SCP1884 
# matrix preprocessed with R function load_Xavier() because it's faster than python ad.read_mtx
# cond is either 'CO' or 'TI'
def preprocess_adata_intestine_Xavier(in_obj, cond):
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
    adata = sc.read_h5ad(in_obj)
    # extract, integrate and rename obs fields
    adata.obs['dataset'] = 'Xavier_' + cond
    adata.obs['tissue'] = 'intestine'
    adata.obs['category'] = 'Unknown'
    df = adata.obs.loc[:,Xavier_fields]
    adata.obs = pandas_rename_columns(df, Xavier_fields, metadata)
    adata.obs['cell_type'] = adata.obs['cell_type'].to_numpy().astype('str')
    
    # rename conditions
    di = {'MONDO_0005011': 'disease', 'PATO_0000461': 'normal'}
    adata.obs['condition'].replace(di, inplace=True)
    
    del adata.raw # this solves the bug in adata.write_h5ad(): "ValueError: '_index' is a reserved name for dataframe columns."
    return adata.copy()


# GSE270788
def preprocess_adata_heart_Amrute(in_dir):
    """
    Put the Henderson adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    Cell ID must be unique 
    
    Parameters
    ----------
    in_dir
        input folder containing the h5 file
    """
    # create adata
    adata = concat_adata_from_multi_h5_heart_Amrute(in_dir, Amrute_healthy)
    adata.obs['condition'] = 'normal'
    adata_cd = concat_adata_from_multi_h5_heart_Amrute(in_dir, Amrute_disease)
    adata_cd.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_cd], merge = "same")
    adata.obs['dataset'] = 'Amrute'
    adata.obs['tissue'] = 'heart'
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['category'] = 'Unknown'
    
    return adata.copy()


# GSE270788
def preprocess_adata_heart_Amrute_filt_by_meta(in_dir):
    """
    Put the Henderson adata format into our format:
    var - ENSEMBL IDs under "gene_ids"
    obs - coarse-grained annotation under "category"
    obs - fine-grained annotation under "cell_type"
    Cell ID must be unique 
    
    Parameters
    ----------
    in_dir
        input folder containing the h5 file
    """

    metadata = in_dir + "/" + Amrute_metadata

    # create adata
    adata = concat_adata_from_multi_h5_heart_Amrute_filt_by_meta(in_dir, Amrute_healthy, metadata)
    adata.obs['condition'] = 'normal'
    adata_cd = concat_adata_from_multi_h5_heart_Amrute_filt_by_meta(in_dir, Amrute_disease, metadata)
    adata_cd.obs['condition'] = 'disease'
    
    # concatenate
    adata = ad.concat([adata, adata_cd], merge = "same")
    adata.obs['dataset'] = 'Amrute'
    adata.obs['tissue'] = 'heart'
    adata.obs['cell_type'] = 'Unknown'
    adata.obs['category'] = 'Unknown'
    
    return adata.copy()


def concat_adata_from_multi_mtx_liver_Henderson(in_dir, prefixes):
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
        sample_name = re.sub(r'GSM[0-9]+_', '', prefix)
        donor_name = re.sub(r'_.+', '', sample_name)
        exp_name = re.sub(r'_.+', '', prefix)
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index
        adata_tmp.obs['library'] = prefix
        adata_tmp.obs['donor_id'] = donor_name
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()


def concat_adata_from_multi_mtx_heart_Du(in_dir, folders):
    """
    Create an anndata object from multiple mtx files
    
    Parameters
    ----------
    in_dir
        input folder containing the matrix folders (*genes*, *barcodes*, *matrix*)
    prefixes
        list of folders, one per matrix
    """
    init = False
    for folder in folders:
        in_folder = in_dir + "/" + folder
        adata_tmp = sc.read_10x_mtx(in_folder, make_unique = True)
        sample_name = re.sub(r'_filtered_feature_bc_matrix', '', folder)
        donor_name = sample_name
        exp_name = re.sub(r'/$', '', folder)
        exp_name = re.sub(r'.+/', '', exp_name)
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index
        adata_tmp.obs['library'] = sample_name
        adata_tmp.obs['donor_id'] = donor_name
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()


def concat_adata_from_multi_mtx_intestine_Sethupathy(in_dir, prefixes):
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
        sample_name = re.sub(r'GSM[0-9]+_', '', prefix)
        donor_name = sample_name
        exp_name = re.sub(r'_.+', '', prefix)
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index
        adata_tmp.obs['library'] = prefix
        adata_tmp.obs['donor_id'] = donor_name
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()


def concat_adata_from_multi_h5_heart_Amrute(in_dir, prefixes):
    """
    Create an anndata object from multiple h5 files
    
    Parameters
    ----------
    in_dir
        input folder containing the files
    prefixes
        list of file prefixes, one per matrix
    """
    init = False
    for prefix in prefixes:
        print(prefix + '\n')
        in_file = in_dir + '/' + prefix + '_filtered_feature_bc_matrix.h5'
        adata_tmp = sc.read_10x_h5(in_file, gex_only = True)
        adata_tmp.var_names_make_unique()
        sample_name = re.sub(r'GSM[0-9]+_', '', prefix)
        donor_name = sample_name
        exp_name = re.sub(r'_.+', '', prefix)
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index # make indexes unique
        adata_tmp.obs['library'] = prefix
        adata_tmp.obs['donor_id'] = donor_name
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()


def concat_adata_from_multi_h5_heart_Amrute_filt_by_meta(in_dir, prefixes, metadata):
    """
    Create an anndata object from multiple h5 files
    
    Parameters
    ----------
    in_dir
        input folder containing the files
    prefixes
        list of file prefixes, one per matrix
    metadata
        input csv file containing metadata information
    """
    df = pd.read_csv(metadata)
    init = False
    for prefix in prefixes:
        in_file = in_dir + '/' + prefix + '_filtered_feature_bc_matrix.h5'
        adata_tmp = sc.read_10x_h5(in_file, gex_only = True)
        adata_tmp.var_names_make_unique()
        sample_name = re.sub(r'GSM[0-9]+_', '', prefix)
        donor_name = sample_name
        exp_name = re.sub(r'_.+', '', prefix)
        
        v = df['sample'] == sample_name
        df_tmp = df.iloc[np.where(v)]
        unnamed_barcodes = [re.sub(r's[0-9]+_', '', i) for i in df_tmp['barcode']] # match barcode names between adata and metadata
        df_tmp.index = unnamed_barcodes # rename indexes so they are the same as in the adata object

        adata_tmp = adata_tmp[unnamed_barcodes,] # subset anndata to match metadata
        adata_tmp.obs = adata_tmp.obs.join(df_tmp) # join based on index
        adata_tmp.obs.index = exp_name + '_' + adata_tmp.obs.index # make indexes unique

        adata_tmp.obs['library'] = prefix
        adata_tmp.obs['donor_id'] = donor_name
        adata_tmp.obs['sample_id'] = sample_name
        if init:
            adata = ad.concat([adata, adata_tmp], merge = "same") # merge = "same" otherwise gene_ids are lost upon concat
        else:
            adata = adata_tmp
            init = True
    return adata.copy()


