
# TODO: save the results of miloR for different parameters in the same object

DIR <- "src/cross_tissue_fibrosis/"

# run DA testing and group adjacent neighbourhood according to FC similarity
runMiloR <- function(q_obj, latent_id, out_miloR_prefix) {
    da_results_file <- paste0(out_miloR_prefix, ".tsv")
    da_results_nhoodgroup_file <- paste0(out_miloR_prefix, "_withNhoodGroups.tsv")
    sce_file <- paste0(out_miloR_prefix, ".Rds")
    cell_nhood_ann <- paste0(out_miloR_prefix, "_nhoodGroup_annotation.tsv")
#    differential_abundance_milo(adata_file = q_obj, latent_id = latent_id, out_prefix = out_miloR_prefix)
    group_nhoods(sce_file = sce_file, da_results_file = da_results_file, out_file = da_results_nhoodgroup_file)
    plot_nhood_groups(da_results_file = da_results_nhoodgroup_file, out_prefix = out_miloR_prefix)
    assign_cells_to_nhood_groups(sce_file = sce_file, da_results_file = da_results_nhoodgroup_file, cell_nhood_ann = cell_nhood_ann)
}

# avg nhood annotation by cell type or category
annMiloR <- function(out_miloR_prefix, q_lab_soft, r_name, meta) {
    da_results_file <- paste0(out_miloR_prefix, ".tsv")
    sce_file <- paste0(out_miloR_prefix, ".Rds")
    da.res <- annotate_neighbourhoods_soft_milo(sce_file = sce_file, da_results_file = da_results_file, soft_ann_file = q_lab_soft)
    da_results_ann_file <- paste0(out_miloR_prefix, "_", r_name, "_", meta, ".tsv")
    write_annotate_neighbourhoods_soft_milo(da.res = da.res, da_results_ann_file = da_results_ann_file)
}

# cell annotation by nhood group 
# marker detection per nhood group
markersMiloR <- function(adata_orig, out_miloR_prefix, meta_files, meta) {
    adata_prefix <- gsub(".h5ad", "", adata_orig)
    out_dea_file <- paste0(out_miloR_prefix, "_nhoodGroup_DEA.tsv")
    out_obj_file <- paste0(out_miloR_prefix, "_SeuratObj.Rds")
    obj <- annotate_object(adata_prefix = adata_prefix, out_miloR_prefix = out_miloR_prefix, meta_files = meta_files)
    for (m in meta) {
#        plot_nhood_group_annotation(meta = obj@meta.data, da_results_file = da_results_nhoodgroup_file, 
#                                    out_prefix = out_miloR_prefix, color.by = m)
        plot_nhood_group_annotation(meta = obj@meta.data, out_miloR_prefix = out_miloR_prefix, color.by = m)
    }
    find_all_nhood_group_markers(obj = obj, out_file = out_dea_file)
    saveRDS(obj, file = out_obj_file)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("\nUsage: <params>\n\n")
    q()
}

PARAMS <- args[1]

### export functions

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
source(file.path(DIR, "differential_abundance.R"))
source(file.path(DIR, "paths.R"))
setwd(WORKING_DIR)

cat("==== Differential abundance ====\n")

### execute

params <- read.table(PARAMS, sep = "=")

folder      <- params[params[,1] == "folder",2]
q_obj_orig  <- params[params[,1] == "q_obj",2] # this obj contains all genes
r_name      <- params[params[,1] == "r_name",2]
q_name      <- params[params[,1] == "q_name",2]
q_batch_key <- params[params[,1] == "q_batch_key",2]
n_top_genes <- params[params[,1] == "n_top_genes",2]
flavor      <- params[params[,1] == "flavor",2]
meta        <- unlist(strsplit(params[params[,1] == "meta",2], split = ","))

q_obj_dir <- create_obj_path(folder, q_name, q_batch_key, n_top_genes, flavor)
q_obj <- create_obj_file(q_obj_dir) # this obj contains only hvgs, selected according to the params above
q_milo_dir <- create_milo_path(folder, q_name, q_batch_key, n_top_genes, flavor)

# scVI
cat("Run miloR on scVI\n")
q_milo <- create_scvi_mod_query(q_milo_dir, r_name)
out_miloR_scvi_prefix <- paste0(q_milo, "_miloR")
latent_id <- paste0("X_scvi_", r_name)
q_lab_hard_1 <- create_labels_file(q_obj_dir, r_name, meta[1])[1]
q_lab_hard_2 <- create_labels_file(q_obj_dir, r_name, meta[2])[1]
meta_files <- c(q_lab_hard_1, q_lab_hard_2)
runMiloR(q_obj = q_obj, latent_id = latent_id, out_miloR_prefix = out_miloR_scvi_prefix)
for (m in meta) {
    # N.B. I can use the cell-type annotation from different embeddings! So do both here
    q_lab <- create_labels_file(q_obj_dir, r_name, m)
    q_lab_soft <- q_lab[2]
    annMiloR(out_miloR_prefix = out_miloR_scvi_prefix, q_lab_soft = q_lab_soft, r_name = r_name, meta = m)
}
markersMiloR(adata_orig = q_obj_orig, out_miloR_prefix = out_miloR_scvi_prefix, meta_files = meta_files, meta = meta)

# scANVI
for (m in meta) {
    cat(paste0("Run miloR on scANVI ", m, "\n"))
    q_milo <- create_scanvi_mod_query(q_milo_dir, r_name, m) 
    out_miloR_scanvi_prefix = paste0(q_milo, "_miloR")
    latent_id <- paste0("X_scanvi_", r_name, "_", m)
    runMiloR(q_obj = q_obj, latent_id = latent_id, out_miloR_prefix = out_miloR_scanvi_prefix)
    q_lab <- create_labels_file(q_obj_dir, r_name, m)
    q_lab_soft <- q_lab[2]
    annMiloR(out_miloR_prefix = out_miloR_scanvi_prefix, q_lab_soft = q_lab_soft, r_name = r_name, meta = m)
    markersMiloR(adata_orig = q_obj_orig, out_miloR_prefix = out_miloR_scanvi_prefix, meta_files = meta_files, meta = meta)
}

q()

