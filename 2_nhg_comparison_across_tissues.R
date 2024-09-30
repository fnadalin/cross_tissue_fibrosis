
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

extract_pair_prefix <- function(params1, params2) {

    folder1      <- params1[params1[,1] == "folder",2]
    r_name1      <- params1[params1[,1] == "r_name",2]
    q_name1      <- params1[params1[,1] == "q_name",2]
    q_batch_key1 <- params1[params1[,1] == "q_batch_key",2]
    n_top_genes1 <- params1[params1[,1] == "n_top_genes",2]
    flavor1      <- params1[params1[,1] == "flavor",2]
    
    folder2      <- params2[params2[,1] == "folder",2]
    r_name2      <- params2[params2[,1] == "r_name",2]
    q_name2      <- params2[params2[,1] == "q_name",2]
    q_batch_key2 <- params2[params2[,1] == "q_batch_key",2]
    n_top_genes2 <- params2[params2[,1] == "n_top_genes",2]
    flavor2      <- params2[params2[,1] == "flavor",2]
    
    if (folder1 != folder2 | q_batch_key1 != q_batch_key2 | n_top_genes1 != n_top_genes2 | flavor1 != flavor2) {
        stop("STOP - Parameters should be consistent\n")
    }

    q_milo_dir <- create_milo_path(folder1, paste(q_name1, q_name2, sep = "-"), q_batch_key1, n_top_genes1, flavor1)
    q_milo <- create_scvi_mod_query(q_milo_dir, paste(r_name1, r_name2, sep = "-"))
    prefix1 <- paste0(q_milo, "_miloR")
    q_milo <- create_scanvi_mod_query(q_milo_dir, paste(r_name1, r_name2, sep = "-"), "category")
    prefix2 <- paste0(q_milo, "_miloR")
    q_milo <- create_scanvi_mod_query(q_milo_dir, paste(r_name1, r_name2, sep = "-"), "cell_type")
    prefix3 <- paste0(q_milo, "_miloR")

    return(c(prefix1,prefix2,prefix3))
} 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <params1> <params2>\n\n")
    q()
}

PARAMS1 <- args[1]
PARAMS2 <- args[2]

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

params1 <- read.table(PARAMS1, sep = "=")
PREFIXES1 <- extract_prefix(params1)

params2 <- read.table(PARAMS2, sep = "=")
PREFIXES2 <- extract_prefix(params2)

PREFIXES_PAIR <- extract_pair_prefix(params1, params2)

for (i in 1:length(PREFIXES1)) {

    PREFIX1 <- PREFIXES1[i]
    PREFIX2 <- PREFIXES2[i]
    PREFIX_PAIR <- PREFIXES_PAIR[i]
    
    nhg1 <- read.table(paste0(PREFIX1, "_withNhoodGroups.tsv"), sep = "\t", header = TRUE)
    nhg2 <- read.table(paste0(PREFIX2, "_withNhoodGroups.tsv"), sep = "\t", header = TRUE)

    q_name1 <- params1[params1[,1] == "q_name",2]
    q_name2 <- params2[params2[,1] == "q_name",2]

    library("reshape2")
    # library("Seurat")
    # library("ggplot2")

    idx <- grep("NhoodGroup", colnames(nhg1))
    if (length(idx) == 0) {
        next
    }
    nhg1[,idx] <- apply(nhg1[,idx,drop=FALSE], 2, as.character)
    idx <- grep("NhoodGroup", colnames(nhg2))
    if (length(idx) == 0) {
        next
    }
    nhg2[,idx] <- apply(nhg2[,idx,drop=FALSE], 2, as.character)

    df1 <- read.table(paste0(PREFIX1, "_nhoodGroup_DEA.tsv"), sep = "\t", header = TRUE)
    markers1 <- read.table(paste0(PREFIX1, "_nhoodGroup_markers.tsv"), sep = "\t", header = TRUE)
    # obj1 <- readRDS(paste0(PREFIX, "_SeuratObj.Rds"))

    df2 <- read.table(paste0(PREFIX2, "_nhoodGroup_DEA.tsv"), sep = "\t", header = TRUE)
    markers2 <- read.table(paste0(PREFIX2, "_nhoodGroup_markers.tsv"), sep = "\t", header = TRUE)
    # obj2 <- readRDS(paste0(PREFIX, "_SeuratObj.Rds"))

    df1$NhoodGroup <- as.character(df1$NhoodGroup)
    df2$NhoodGroup <- as.character(df2$NhoodGroup)

    nhg_params1 <- unique(df1$NhoodGroupParams)
    nhg_params2 <- unique(df2$NhoodGroupParams)

    row_names <- unlist(lapply(nhg_params1, function(x) paste(x, unique(df1$NhoodGroup[df1$NhoodGroupParams == x]), sep = ".")))
    col_names <- unlist(lapply(nhg_params2, function(x) paste(x, unique(df2$NhoodGroup[df2$NhoodGroupParams == x]), sep = ".")))
    jaccard <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
    rownames(jaccard) <- row_names
    colnames(jaccard) <- col_names
    shared_markers <- mean_auc <- mean_logFC <- delta_logFC <- jaccard

    count1 <- 0
    for (param1 in nhg_params1) {

        # compute markers
        gene.list.1 <- lapply(markers1$group[markers1$param == param1], function(x) 
                          unlist(strsplit(markers1$markers[markers1$group == x], split = ",")))
        names(gene.list.1) <- as.character(markers1$group[markers1$param == param1])
        n <- length(gene.list.1)

        I <- count1 + 1:n
        i_labels <- unique(df1$NhoodGroup[df1$NhoodGroupParams == param1])
        
        count2 <- 0
        for (param2 in nhg_params2) {

        # compute markers
        gene.list.2 <- lapply(markers2$group[markers2$param == param2], function(x) 
                              unlist(strsplit(markers2$markers[markers2$group == x], split = ",")))
        names(gene.list.2) <- as.character(markers2$group[markers2$param == param2])
        m <- length(gene.list.2)
        
        J <- count2 + 1:m
        j_labels <- unique(df2$NhoodGroup[df2$NhoodGroupParams == param2])

        # compute Jaccard between nhoodGroups
        M <- Jaccard(gene.list.1 = gene.list.1, gene.list.2 = gene.list.2)
        jaccard[I,J] <- M[i_labels, j_labels] # N.B.: jaccard is usually very low, AUC is very high, so maybe not very useful to compare 
                                              # perhaps better to normalize
                                              
        shared_markers[I,J] <- shared_marker_matrix(gene.list.1, gene.list.2)
        
        # compute markers AUC
#        M <- mean_shared_markers_auc_parallel(obj1 = obj1, obj2 = obj2, param1 = param1, param2 = param2, 
#                                              gene.list.1 = gene.list.1, gene.list.2 = gene.list.2)
#        mean_auc[I,J] <- M[i_labels, j_labels]
        
        # compute the average logFC
        M <- mean_abs_logFC_nhood_groups_from_data_frame(df1 = markers1, df2 = markers2, param1 = param1, param2 = param2)
        mean_logFC[I,J] <- M[i_labels, j_labels] 
        
        # compute the logFC similarity
        M <- delta_logFC_nhood_groups_from_data_frame(df1 = markers1, df2 = markers2, param1 = param1, param2 = param2)
        delta_logFC[I,J] <- M[i_labels, j_labels] 

        count2 <- count2 + m
        }
        count1 <- count1 + n
    }

    # normalise
    mean_logFC_norm <- (mean_logFC-min(mean_logFC))/(max(mean_logFC)-min(mean_logFC))
    logFC_sim_norm <- 1 - (delta_logFC-min(delta_logFC))/(max(delta_logFC)-min(delta_logFC))
    jaccard_norm <- (jaccard-min(jaccard))/(max(jaccard)-min(jaccard))
    colnames(mean_logFC_norm) <- colnames(logFC_sim_norm) <- colnames(jaccard_norm) <- col_names
    rownames(mean_logFC_norm) <- rownames(logFC_sim_norm) <- rownames(jaccard_norm) <- row_names

    # combine the scores
    F1_logFC <- F1_score(x = mean_logFC_norm, y = logFC_sim_norm)
    colnames(F1_logFC) <- col_names
    rownames(F1_logFC) <- row_names
    combined1 <- F1_score(jaccard_norm, F1_logFC)
    combined2 <- F1_score(jaccard_norm, logFC_sim_norm)
    combined3 <- (jaccard_norm + F1_logFC)/2
    combined4 <- (jaccard_norm + logFC_sim_norm)/2
    colnames(combined1) <- colnames(combined2) <- colnames(combined3) <- colnames(combined4) <- col_names
    rownames(combined1) <- rownames(combined2) <- rownames(combined3) <- rownames(combined4) <- row_names

    # print the matrices
    write.table(jaccard, file = paste0(PREFIX_PAIR, "_jaccard.tsv"), sep = "\t", quote = FALSE)
    write.table(mean_logFC, file = paste0(PREFIX_PAIR, "_mean_logFC.tsv"), sep = "\t", quote = FALSE)
    write.table(logFC_sim_norm, file = paste0(PREFIX_PAIR, "_logFC_sim.tsv"), sep = "\t", quote = FALSE)
    write.table(F1_logFC, file = paste0(PREFIX_PAIR, "_F1_logFC.tsv"), sep = "\t", quote = FALSE)
    write.table(shared_markers, file = paste0(PREFIX_PAIR, "_shared_markers.tsv"), sep = "\t", quote = FALSE)
    write.table(combined1, file = paste0(PREFIX_PAIR, "_combined1_score.tsv"), sep = "\t", quote = FALSE)
    write.table(combined2, file = paste0(PREFIX_PAIR, "_combined2_score.tsv"), sep = "\t", quote = FALSE)
    write.table(combined3, file = paste0(PREFIX_PAIR, "_combined3_score.tsv"), sep = "\t", quote = FALSE)
    write.table(combined4, file = paste0(PREFIX_PAIR, "_combined4_score.tsv"), sep = "\t", quote = FALSE)
}

q()


