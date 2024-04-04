
create_obj_path <- function(parent_dir, name, batch_key, n_top_genes, flavor) {
    dirname <- file.path(parent_dir, "adata", name, batch_key, paste0(flavor, "_", n_top_genes))
    dir.create(dirname, recursive = TRUE, showWarnings = FALSE) 
    return(dirname)
}


create_obj_file <- function(obj_dir) {
    return(file.path(obj_dir, "adata.h5ad"))
}


create_milo_path <- function(parent_dir, name, batch_key, n_top_genes = 5000, flavor = "seurat_v3") {
    dirname <- file.path(parent_dir, "differential_analysis", name, batch_key, paste0(flavor, "_", n_top_genes))
    dir.create(dirname, recursive = TRUE, showWarnings = FALSE) 
    return(dirname)
}


create_scvi_mod_query <- function(q_mod_dir, ref_name) {
    return(file.path(q_mod_dir, paste0("scvi_", ref_name)))
}


create_scanvi_mod_query <- function(q_mod_dir, ref_name, meta) {
    return(file.path(q_mod_dir, paste0("scanvi_", ref_name, "_", meta)))
}


create_labels_file <- function(q_obj_dir, ref_name, meta) {
    q_lab_hard <- file.path(q_obj_dir, paste0("scanvi_", ref_name, "_", meta, "_hard_labels.tsv"))
    q_lab_soft <- file.path(q_obj_dir, paste0("scanvi_", ref_name, "_", meta, "_soft_labels.tsv"))
    return(c(q_lab_hard, q_lab_soft))
}


