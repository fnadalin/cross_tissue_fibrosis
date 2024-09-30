
# only extract file names

DIR <- "src/cross_tissue_fibrosis/"

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

    q_obj_dir <- create_obj_path(folder, q_name, q_batch_key, n_top_genes, flavor)
    q_lab_cat <- create_labels_file(q_obj_dir, r_name, "category")[1]
    q_lab_ct <- create_labels_file(q_obj_dir, r_name, "cell_type")[1]

    return(c(prefix1,prefix2,prefix3,q_lab_cat,q_lab_ct))
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
if (length(args) < 1) {
    cat("\nUsage: <out_prefix> <params1> [... <paramsN>]\n\n")
    q()
}

OUT_PREFIX <- args[1]
PARAMS <- args[2:length(args)]

### export functions

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
source(file.path(DIR, "paths.R"))
setwd(WORKING_DIR)

cat("==== File generation ====\n")

### generate annotation file table

df_ann_1 <- as.data.frame(matrix(NA, nrow = length(PARAMS), ncol = 4))
colnames(df_ann_1) <- c("tissue.id","cell.nhg","cell.category","cell.cell_type")
df_ann_2 <- df_ann_3 <- df_ann_1
# NEW: consider only the tissues with at least one nhg found
PARAM_ID_SCVI <- c()
PARAM_ID_SCANVI_CAT <- c()
PARAM_ID_SCANVI_CT <- c()
for (i in 1:length(PARAMS)) {

    params <- read.table(PARAMS[i], sep = "=")
    PREFIXES <- extract_prefix(params)
    
    field1 <- params[params[,1] == "q_name",2]
    field3 <- PREFIXES[4]
    field4 <- PREFIXES[5]
    
    if (file.exists(paste0(PREFIXES[1], "_nhoodGroup_markers.tsv")) {
        field2 <- paste0(PREFIXES[1], "_nhoodGroup_annotation.tsv")
        df_ann_1[i,] <- c(field1, field2, field3, field4)
        PARAM_ID_SCVI <- c(PARAM_ID_SCVI, i)
    }
    if (file.exists(paste0(PREFIXES[2], "_nhoodGroup_markers.tsv")) {
        field2 <- paste0(PREFIXES[2], "_nhoodGroup_annotation.tsv")
        df_ann_2[i,] <- c(field1, field2, field3, field4)
        PARAM_ID_SCANVI_CAT <- c(PARAM_ID_SCANVI_CAT, i)
    }
    if (file.exists(paste0(PREFIXES[3], "_nhoodGroup_markers.tsv")) {
        field2 <- paste0(PREFIXES[3], "_nhoodGroup_annotation.tsv")
        df_ann_3[i,] <- c(field1, field2, field3, field4)
        PARAM_ID_SCANVI_CT <- c(PARAM_ID_SCANVI_CT, i)
    }
}
write.table(df_ann_1, file = paste0(OUT_PREFIX, "_scvi_annotation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_ann_2, file = paste0(OUT_PREFIX, "_scanvi_category_annotation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_ann_3, file = paste0(OUT_PREFIX, "_scanvi_cell_type_annotation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

### generate pairwise comparison file table

for (c in 1:4) {

    df_ann_1 <- as.data.frame(matrix(NA, nrow = length(PARAM_ID_SCVI)*(length(PARAM_ID_SCVI)-1)/2, ncol = 5))
    df_ann_2 <- as.data.frame(matrix(NA, nrow = length(PARAM_ID_SCANVI_CAT)*(length(PARAM_ID_SCANVI_CAT)-1)/2, ncol = 5))
    df_ann_3 <- as.data.frame(matrix(NA, nrow = length(PARAM_ID_SCANVI_CT)*(length(PARAM_ID_SCANVI_CT)-1)/2, ncol = 5))
    colnames(df_ann_1) <- colnames(df_ann_2) <- colnames(df_ann_3) <- c("tissue1.id","tissue2.id","nhg.1","nhg.2","pairwise.score")
    
    k <- 1
    for (i in PARAM_ID_SCVI[1:(length(PARAM_ID_SCVI)-1))) {
    
        params1 <- read.table(PARAMS[i], sep = "=")
        PREFIXES1 <- extract_prefix(params1)
        
        field1 <- params1[params1[,1] == "q_name",2]
        field3 <- paste0(PREFIXES1[1], "_nhoodGroup_markers.tsv")
        
        for (j in PARAM_ID_SCVI[(i+1):length(PARAM_ID_SCVI)]) {
    
            params2 <- read.table(PARAMS[j], sep = "=")
            PREFIXES2 <- extract_prefix(params2)
            PAIR_PREFIXES <- extract_pair_prefix(params1, params2)
        
            field2 <- params2[params2[,1] == "q_name",2]
            field4 <- paste0(PREFIXES2[1], "_nhoodGroup_markers.tsv")
            field5 <- paste0(PAIR_PREFIXES[1], "_combined", c, "_score.tsv")
            df_ann_1[k,] <- c(field1, field2, field3, field4, field5)
            
            k <- k+1
        }
    }
    write.table(df_ann_1, file = paste0(OUT_PREFIX, "_scvi_combined", c, "_score.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

    k <- 1
    for (i in 1:(length(PARAMS)-1)) {
    
        params1 <- read.table(PARAMS[i], sep = "=")
        PREFIXES1 <- extract_prefix(params1)
        
        field1 <- params1[params1[,1] == "q_name",2]
        field3 <- paste0(PREFIXES1[2], "_nhoodGroup_markers.tsv")
        
        for (j in (i+1):length(PARAMS)) {
    
            params2 <- read.table(PARAMS[j], sep = "=")
            PREFIXES2 <- extract_prefix(params2)
            PAIR_PREFIXES <- extract_pair_prefix(params1, params2)
        
            field2 <- params2[params2[,1] == "q_name",2]
            field4 <- paste0(PREFIXES2[2], "_nhoodGroup_markers.tsv")
            field5 <- paste0(PAIR_PREFIXES[2], "_combined", c, "_score.tsv")
            df_ann_2[k,] <- c(field1, field2, field3, field4, field5)
            
            k <- k+1
        }
    }    
    write.table(df_ann_2, file = paste0(OUT_PREFIX, "_scanvi_category_combined", c, "_score.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

    k <- 1
    for (i in 1:(length(PARAMS)-1)) {
    
        params1 <- read.table(PARAMS[i], sep = "=")
        PREFIXES1 <- extract_prefix(params1)
        
        field1 <- params1[params1[,1] == "q_name",2]
        field3 <- paste0(PREFIXES1[3], "_nhoodGroup_markers.tsv")
        
        for (j in (i+1):length(PARAMS)) {
    
            params2 <- read.table(PARAMS[j], sep = "=")
            PREFIXES2 <- extract_prefix(params2)
            PAIR_PREFIXES <- extract_pair_prefix(params1, params2)
        
            field2 <- params2[params2[,1] == "q_name",2]
            field4 <- paste0(PREFIXES2[3], "_nhoodGroup_markers.tsv")
            field5 <- paste0(PAIR_PREFIXES[3], "_combined", c, "_score.tsv")
            df_ann_3[k,] <- c(field1, field2, field3, field4, field5)
            
            k <- k+1
        }
    }  
    write.table(df_ann_3, file = paste0(OUT_PREFIX, "_scanvi_cell_type_combined", c, "_score.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

q()


