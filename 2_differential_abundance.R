
DIR <- "src/cross_tissue_fibrosis/"

runMiloR <- function(q_obj, latent_id, out_miloR_prefix) {
    da_results_file <- paste0(out_miloR_prefix, ".tsv")
    differential_abundance_milo(q_obj, latent_id = latent_id, out_miloR_prefix) 
}

annMiloR <- function(out_miloR_prefix, q_lab_soft, r_name, meta) {
    sce_file <- paste0(out_miloR_prefix, ".Rds")
    da_results_file <- paste0(out_miloR_prefix, ".tsv")
    da.res <- annotate_neighbourhoods_soft_milo(sce_file, da_results_file, q_lab_soft)
    da_results_ann_file <- paste0(out_miloR_prefix, "_", r_name, "_", meta, ".tsv")
    write_annotate_neighbourhoods_soft_milo(da.res, da_results_ann_file)
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

### execute

params <- read.table(PARAMS, sep = "=")

folder      <- params[params[,1] == "folder",2]
q_obj       <- params[params[,1] == "q_obj",2]
r_name      <- params[params[,1] == "r_name",2]
q_name      <- params[params[,1] == "q_name",2]
q_batch_key <- params[params[,1] == "q_batch_key",2]
n_top_genes <- params[params[,1] == "n_top_genes",2]
flavor      <- params[params[,1] == "flavor",2]
meta        <- unlist(strsplit(params[params[,1] == "meta",2], split = ","))

q_obj_dir <- create_obj_path(folder, q_name, q_batch_key, n_top_genes, flavor)
q_obj <- create_obj_file(q_obj_dir)
q_milo_dir <- create_milo_path(folder, q_name, q_batch_key, n_top_genes, flavor)

# scVI
q_milo <- create_scvi_mod_query(q_milo_dir, r_name)
out_miloR_scvi_prefix <- paste0(q_milo, "_miloR")
latent_id <- paste0("X_scvi_", r_name)
runMiloR(q_obj, latent_id = latent_id, out_miloR_scvi_prefix)
for (m in meta) {
    # N.B. I can use the cell-type annotation from different embeddings! So do both here
    q_lab <- create_labels_file(q_obj_dir, r_name, m)
    q_lab_soft <- q_lab[2]
    annMiloR(out_miloR_scvi_prefix, q_lab_soft, r_name, m)
}

# scANVI
for (m in meta) {
    q_milo <- create_scanvi_mod_query(q_milo_dir, r_name, m) 
    out_miloR_scanvi_prefix = paste0(q_milo, "_miloR")
    latent_id <- paste0("X_scanvi_", r_name, "_", m)
    runMiloR(q_obj, latent_id = latent_id, out_miloR_scanvi_prefix)
    q_lab <- create_labels_file(q_obj_dir, r_name, m)
    q_lab_soft <- q_lab[2]
    annMiloR(out_miloR_scanvi_prefix, q_lab_soft, r_name, m)
}

q()

