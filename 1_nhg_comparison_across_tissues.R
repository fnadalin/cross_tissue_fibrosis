
DIR <- "src/cross_tissue_fibrosis/"
LOGFC_MARKER <- 1

extract_prefix <- function(params) {

    folder      <- params[params[,1] == "folder",2]
    r_name      <- params[params[,1] == "r_name",2]
    q_name      <- params[params[,1] == "q_name",2]
    q_batch_key <- params[params[,1] == "q_batch_key",2]
    n_top_genes <- params[params[,1] == "n_top_genes",2]
    flavor      <- params[params[,1] == "flavor",2]

    q_milo_dir <- create_milo_path(folder, q_name, q_batch_key, n_top_genes, flavor)
    q_milo <- create_scvi_mod_query(q_milo_dir, r_name)
    prefix1 <- paste0(q_milo, "_miloR")
    q_milo <- create_scanvi_mod_query(q_milo_dir, r_name, "category")
    prefix2 <- paste0(q_milo, "_miloR")
    q_milo <- create_scanvi_mod_query(q_milo_dir, r_name, "cell_type")
    prefix3 <- paste0(q_milo, "_miloR")

    return(c(prefix1,prefix2,prefix3))
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
source(file.path(DIR, "paths.R"))
source(file.path(DIR, "cross_tissue_comparison.R"))
setwd(WORKING_DIR)

cat("==== Nhood comparison ====\n")

### execute

params <- read.table(PARAMS, sep = "=")
PREFIXES <- extract_prefix(params)

for (PREFIX in PREFIXES) {

    df <- read.table(paste0(PREFIX, "_nhoodGroup_DEA.tsv"), sep = "\t", header = TRUE)
    nhg <- read.table(paste0(PREFIX, "_withNhoodGroups.tsv"), sep = "\t", header = TRUE)
    # obj <- readRDS(paste0(PREFIX, "_SeuratObj.Rds"))

    q_name <- params[params[,1] == "q_name",2]

    library("reshape2")
    # library("Seurat")
    # library("ggplot2")

    df$NhoodGroup <- as.character(df$NhoodGroup)

    idx <- grep("hood", colnames(nhg))
    nhg[,idx] <- apply(nhg[,idx], 2, as.character)

    nhg_params <- unique(df$NhoodGroupParams)
    param_names <- unlist(lapply(nhg_params, function(x) rep(x, length(unique(df$NhoodGroup[df$NhoodGroupParams == x])))))
    group_names <- unlist(lapply(nhg_params, function(x) unique(df$NhoodGroup[df$NhoodGroupParams == x])))

    markers_all <- c()
    nhg_logfc_all <- c()

    # N.B. consider only the nhg that have some markers associated with them
    for (param in nhg_params) {

        i_labels <- unique(df$NhoodGroup[df$NhoodGroupParams == param])
        
        l <- extract_nhood_group_markers(df = df, param = param, LOGFC_MARKER = LOGFC_MARKER)
        markers <- unlist(lapply(l, function(x) paste(x, collapse = ",")))
        
        nhg_logfc <- logFC_nhood_groups(nhg = nhg, param = param)
        
        markers_all <- c(markers_all, markers[i_labels])
        nhg_logfc_all <- c(nhg_logfc_all, nhg_logfc[i_labels])
    }

    df_nhg <- data.frame(param = param_names, group = group_names, logfc = nhg_logfc_all, markers = markers_all)
    write.table(df_nhg, file = paste0(PREFIX, "_nhoodGroup_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}


q()


