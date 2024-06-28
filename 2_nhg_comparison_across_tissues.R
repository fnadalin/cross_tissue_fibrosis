
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
    prefix <- paste0(q_milo, "_miloR")

    return(prefix)
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
    prefix <- paste0(q_milo, "_miloR")

    return(prefix)
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
PREFIX <- extract_prefix(params1)
df1 <- read.table(paste0(PREFIX, "_nhoodGroup_DEA.tsv"), sep = "\t", header = TRUE)
nhg1 <- read.table(paste0(PREFIX, "_withNhoodGroups.tsv"), sep = "\t", header = TRUE)
markers1 <- read.table(paste0(PREFIX, "_nhoodGroup_markers.tsv"), sep = "\t", header = TRUE)
# obj1 <- readRDS(paste0(PREFIX, "_SeuratObj.Rds"))

params2 <- read.table(PARAMS2, sep = "=")
PREFIX <- extract_prefix(params2)
df2 <- read.table(paste0(PREFIX, "_nhoodGroup_DEA.tsv"), sep = "\t", header = TRUE)
nhg2 <- read.table(paste0(PREFIX, "_withNhoodGroups.tsv"), sep = "\t", header = TRUE)
markers2 <- read.table(paste0(PREFIX, "_nhoodGroup_markers.tsv"), sep = "\t", header = TRUE)
# obj2 <- readRDS(paste0(PREFIX, "_SeuratObj.Rds"))

q_name1 <- params1[params1[,1] == "q_name",2]
q_name2 <- params2[params2[,1] == "q_name",2]

PREFIX_PAIR <- extract_pair_prefix(params1, params2)

library("reshape2")
# library("Seurat")
# library("ggplot2")

df1$NhoodGroup <- as.character(df1$NhoodGroup)
df2$NhoodGroup <- as.character(df2$NhoodGroup)

idx <- grep("hood", colnames(nhg1))
nhg1[,idx] <- apply(nhg1[,idx], 2, as.character)
idx <- grep("hood", colnames(nhg2))
nhg2[,idx] <- apply(nhg2[,idx], 2, as.character)

params1 <- unique(df1$NhoodGroupParams)
params2 <- unique(df2$NhoodGroupParams)

row_names <- unlist(lapply(params1, function(x) paste(x, unique(df1$NhoodGroup[df1$NhoodGroupParams == x]), sep = ".")))
col_names <- unlist(lapply(params2, function(x) paste(x, unique(df2$NhoodGroup[df2$NhoodGroupParams == x]), sep = ".")))
jaccard <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
rownames(jaccard) <- row_names
colnames(jaccard) <- col_names
combined <- shared_markers <- mean_auc <- mean_logFC <- logFC_sim <- F1_markers <- F1_logFC <- jaccard

count1 <- 0
for (param1 in params1) {

    # compute markers
    gene.list.1 <- lapply(markers1$group[markers1$param == param1], function(x) 
                          unlist(strsplit(markers1$markers[markers1$group == x], split = ",")))
    names(gene.list.1) <- as.character(markers1$group[markers1$param == param1])
    n <- length(gene.list.1)

    I <- count1 + 1:n
    i_labels <- unique(df1$NhoodGroup[df1$NhoodGroupParams == param1])
    
    count2 <- 0
    for (param2 in params2) {

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
        mean_logFC[I,J] <- M[i_labels, j_labels] # TODO: normalise wrt to ALL the elements of the matrix
        
        # compute the logFC similarity
        M <- logFC_similarity_nhood_groups_from_data_frame(df1 = markers1, df2 = markers2, param1 = param1, param2 = param2)
        logFC_sim[I,J] <- M[i_labels, j_labels] # TODO: normalise wrt to ALL elements of the matris
        
        # calculate F1 scores
#        F1_markers[I,J] <- F1_score(x = jaccard[I,J], y = mean_auc[I,J])
        F1_logFC[I,J] <- F1_score(x = mean_logFC[I,J], y = logFC_sim[I,J]) # TODO: put normalisation and F1 calculation in a different for loop 
        
        count2 <- count2 + m
    }
    count1 <- count1 + n
}

# combine the scores
jaccard_norm <- (jaccard-min(jaccard))/(max(jaccard)-min(jaccard))
F1_logFC_norm <- (F1_logFC-min(F1_logFC))/(max(F1_logFC)-min(F1_logFC))
combined <- F1_score(jaccard_norm, F1_logFC_norm)

# print the matrices
write.table(jaccard, file = paste0(PREFIX_PAIR, "_jaccard.tsv"), sep = "\t", quote = FALSE)
write.table(mean_logFC, file = paste0(PREFIX_PAIR, "_mean_logFC.tsv"), sep = "\t", quote = FALSE)
write.table(logFC_sim, file = paste0(PREFIX_PAIR, "_logFC_sim.tsv"), sep = "\t", quote = FALSE)
write.table(F1_logFC, file = paste0(PREFIX_PAIR, "_F1_logFC.tsv"), sep = "\t", quote = FALSE)
write.table(combined, file = paste0(PREFIX_PAIR, "_combined_score.tsv"), sep = "\t", quote = FALSE)
write.table(shared_markers, file = paste0(PREFIX_PAIR, "_shared_markers.tsv"), sep = "\t", quote = FALSE)


q()


