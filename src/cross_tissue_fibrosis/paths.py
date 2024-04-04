### environment

import os

### functions

def create_scvi_mod_ref(ref_mod_dir):
    return os.path.join(ref_mod_dir, "scvi")


def create_scvi_mod_query(q_mod_dir, ref_name):
    return os.path.join(q_mod_dir, "scvi_" + ref_name)


def create_scanvi_mod_ref(ref_mod_dir, meta):
    return os.path.join(ref_mod_dir, "scanvi_" + meta)


def create_scanvi_mod_query(q_mod_dir, ref_name, meta):
    return os.path.join(q_mod_dir, "scanvi_" + ref_name + "_" + meta)


def create_obj_file(obj_dir):
    return os.path.join(obj_dir, "adata.h5ad")


def create_labels_file(q_obj_dir, ref_name, meta):
    q_lab_hard = os.path.join(q_obj_dir, "scanvi_" + ref_name + "_" + meta + "_hard_labels.tsv")
    q_lab_soft = os.path.join(q_obj_dir, "scanvi_" + ref_name + "_" + meta + "_soft_labels.tsv")
    return q_lab_hard, q_lab_soft


def create_mod_path(parent_dir, name, batch_key, n_top_genes = 5000, flavor = "seurat_v3"):       
    dirname = os.path.join(parent_dir, "models", name, batch_key, flavor + "_" + str(n_top_genes))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)    
    return dirname


def create_obj_path(parent_dir, name, batch_key, n_top_genes = 5000, flavor = "seurat_v3"):    
    dirname = os.path.join(parent_dir, "adata", name, batch_key, flavor + "_" + str(n_top_genes))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)  
    return dirname


def create_milo_path(parent_dir, name, batch_key, n_top_genes = 5000, flavor = "seurat_v3"):
    dirname = os.path.join(parent_dir, "differential_analysis", name, batch_key, flavor + "_" + str(n_top_genes))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)  
    return dirname




