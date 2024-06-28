
# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

library("anndataR")
library("SingleCellExperiment")
library("miloR")
library("dplyr")
library("scater") # only for logNormCounts()
library("ggplot2")
library("Matrix")
library("ggplot2")
library("Seurat")
library("SeuratDisk")

# MAX_LOGFC_DELTA <- 0.5
# NH_OVERLAP <- 10
# MIN_PERC <- 0.7
# TEST_USE <- "MAST"

set.seed(123)

# FOR TESTING: use the atlas sample dataset with "disease" assigned to the last 500 cells: tests/atlas_liver_sample_simConditions.h5ad
# (in the query sample dataset, no DA nhood was found)
# N.B.: make sure that cells in the same sample are all assigned to the same condition!!!
differential_abundance_milo <- function(adata_file, latent_id = 'X_scvi', out_prefix) {
    
    # convert adata to SingleCellExperiment
    adata <- read_h5ad(adata_file, to = "InMemoryAnnData")
    sce <- adata$to_SingleCellExperiment()
    
    # prepare the SingleCellExperiment object for miloR
    colnames(sce) <- adata$obs_names # NEW: assign cell names
    names(assays(sce)) <- "counts"
    reducedDim(sce) <- adata$obsm[[latent_id]]
    reducedDimNames(sce) <- latent_id
    d <- ncol(reducedDim(sce))
    
    sce <- Milo(sce)
    n_controls <- length(unique(adata$obs[["sample_id"]][adata$obs[["condition"]] == "normal"]))
    n_querys <- length(unique(adata$obs[["sample_id"]][adata$obs[["condition"]] == "disease"]))
    #  Set max to 200 or memory explodes for large datasets
    k = min((n_controls + n_querys) * 5, 200)
    sce <- buildGraph(sce, k = k, d = d, reduced.dim = latent_id)
    sce <- makeNhoods(sce, prop = 0.1, k = k, d = d, reduced_dims = latent_id, refined = TRUE)
    
    sce <- countCells(sce, meta.data = as.data.frame(colData(sce)), sample="sample_id")
    
    design <- data.frame(colData(sce))[,c("sample_id", "condition")]
    design <- distinct(design)
    rownames(design) <- design$sample_id
    ## Reorder rownames to match columns of nhoodCounts(milo)
    design <- design[colnames(nhoodCounts(sce)), , drop=FALSE]
    
    sce <- calcNhoodDistance(sce, reduced.dim = latent_id, d = d)
    da_results <- testNhoods(sce, design = ~condition, design.df = design, reduced.dim = latent_id)
    
    out_da <- paste0(out_prefix, ".tsv")
    write.table(da_results, file = out_da, quote = FALSE, sep = "\t", row.names = FALSE)
    
    out_sce <- paste0(out_prefix, ".Rds")
    saveRDS(sce, file = out_sce) 
}


# adapted from annotateNhoods from package:miloR
# compute a nhood label probability by majority voting
annotate_neighbourhoods_soft_milo <- function(sce_file, da_results_file, soft_ann_file) {

    x <- readRDS(sce_file)
    da.res <- read.table(da_results_file, sep = "\t", header = TRUE)
    soft_ann <- read.table(soft_ann_file, sep = "\t", header = TRUE, row.names = 1)

    if (!is(x, "Milo")) {
        stop("Unrecognised input type - must be of class Milo")
    }
    if (ncol(nhoods(x)) != nrow(da.res)) {
        stop("the number of rows in da.res does not match the number of neighbourhoods in nhoods(x). Are you sure da.res is the output of testNhoods(x)?")
    }
    if (nrow(nhoods(x)) != nrow(soft_ann)) {
        stop("the number of rows in soft_ann does not match the number of cells in nhoods(x). Are you sure Milo was run on the same object?")
    }

    nhood_counts <- vapply(seq_len(ncol(nhoods(x))), FUN = function(n) colSums(soft_ann[which(nhoods(x)[,n] == 1),]), FUN.VALUE = numeric(ncol(soft_ann)))
    nhood_counts <- nhood_counts/colSums(nhood_counts)
    nhood_counts <- t(nhood_counts)
    rownames(nhood_counts) <- seq_len(ncol(nhoods(x)))
    colnames(nhood_counts) <- colnames(soft_ann)
    da.res <- cbind(da.res, nhood_counts)
    
    return(da.res)
}


# output the table with da.res + annotation ready for python
write_annotate_neighbourhoods_soft_milo <- function(da.res, da_results_ann_file) {

    df <- data.frame(rownames(da.res))
    da.res <- cbind(df, da.res)
    colnames(da.res)[1] <- ""
    write.table(da.res, file = da_results_ann_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


group_nhoods <- function(sce_file, da_results_file, out_file, MAX_LOGFC_DELTA = c(0.5,1,2), NH_OVERLAP = c(1,10,20,50), FDR = 0.1) {

    # TODO: debug the code below - OOM for heart and lung!!!
    sce <- readRDS(sce_file)
    da.res <- read.table(da_results_file, sep = "\t", header = TRUE)
    idx <- which(da.res$SpatialFDR <= FDR) # NEW: only consider significantly DA nhoods
    for (max_logfc_delta in MAX_LOGFC_DELTA) {
        for (nh_overlap in NH_OVERLAP) {
#            da.res <- groupNhoods(sce, da.res, max.lfc.delta = max_logfc_delta, overlap = nh_overlap, 
#                                  da.fdr = 1, merge.discord = FALSE)
            # N.B. da.fdr does not select nhoohd for grouping, but sets a criterion on which nhoods should be considered for logFC consistency
            da.res <- groupNhoods(sce, da.res, max.lfc.delta = max_logfc_delta, overlap = nh_overlap, 
                                  subset.nhoods = idx, da.fdr = 1, merge.discord = FALSE, compute.new = TRUE)
            colnames(da.res)[colnames(da.res) == "NhoodGroup"] <- paste0("NhoodGroup_logFC", max_logfc_delta, "_overlap", nh_overlap)
        }
    }
    write.table(da.res, file = out_file, sep = "\t", quote = FALSE)
}


plot_nhood_groups <- function(da_results_file, out_prefix) {

    da.res <- read.table(da_results_file, sep = "\t", header = TRUE)
    fields <- grep("NhoodGroup", colnames(da.res), value = TRUE)
    for (field in fields) {
        df <- da.res
        colnames(df)[colnames(df) == field] <- "NhoodGroup"
        df <- df[!is.na(df$NhoodGroup),] # remove the nhood group formed by NA nhoods (the ones discarded in group_nhoods())
        groups <- unique(df$NhoodGroup)
        avg_logfc <- unlist(lapply(groups, function(x) mean(df$logFC[df$NhoodGroup == x])))
        groups <- groups[order(avg_logfc)]
        df$NhoodGroup <- factor(df$NhoodGroup, levels = groups)

        # plot nhood groups
        pdf(paste0(out_prefix, "_", field, "_beeswarmDA.pdf"), height = 10)
        print(plotDAbeeswarm(df, group.by = "NhoodGroup") + ggtitle(field))
        dev.off()
    }
}


assign_cells_to_nhood_groups <- function(sce_file, da_results_file, cell_nhood_ann, MIN_PERC = 0.7) {

    sce <- readRDS(sce_file)
    da.res <- read.table(da_results_file, sep = "\t", header = TRUE)
    
    meta <- as.data.frame(matrix(NA, nrow = nrow(sce@nhoods), ncol = 0))
    rownames(meta) <- colnames(sce)
    fields <- grep("NhoodGroup", colnames(da.res), value = TRUE)
    for (field in fields) {
        df <- da.res
        colnames(df)[colnames(df) == field] <- "NhoodGroup"
        
        # create cellxnhoodgroup matrix
        # sce@nhoods: rows = cells, columns = nhoods
        # df$X: nhood ID
        # df$NhoodGroup: nhood group
        nhood_idx <- which(!is.na(df$NhoodGroup))
        df <- df[nhood_idx,]
        groups <- unique(df$NhoodGroup)
        M1 <- sapply(groups, function(x) as.numeric(df$NhoodGroup == x))
        M1 <- Matrix(M1) # nhoodxnhoodgroup matrix
        rownames(M1) <- rownames(df)
        colnames(M1) <- groups
#        M <- sce@nhoods %*% M1 # cellxnhoodgroup matrix containing the count of nhood
        M <- sce@nhoods[,nhood_idx] %*% M1 # NEW: subset the sce to only contain non-NA nhoods
    
        # assign each cell to a group by selecting the one with highest number of nhood containing the cell
        # delete the assignment if the top groups contain less than MIN_PERC fraction of nhoods that contain the cell
        assignment <- apply(M, 1, which.max)
        max.val <- apply(M/rowSums(M), 1, max)
        idx <- which(max.val < MIN_PERC)
        assignment[idx] <- rep(NaN, length(idx))
        
        cnames <- colnames(meta)
        meta <- cbind(meta, x = assignment)
        colnames(meta) <- c(cnames, field)
    }
    write.table(meta, file = cell_nhood_ann, sep = "\t", quote = FALSE)
}


# the object should contain *all* genes 
# the meta_file_list is a vector of file names that contain metadata for the count matrix
# ideally, it should be: nhood_group metadata (cell assignment to nhood_groups), hard cell-type annotation, hard category annotation
annotate_object_OLD <- function(adata_prefix, out_miloR_prefix, meta_files = NULL) {

    adata_file <- paste0(adata_prefix, ".h5ad")
    seurat_file <- paste0(adata_prefix, ".h5seurat")
    adata <- read_h5ad(adata_file, to = "InMemoryAnnData")
    Convert(adata_file, dest = "h5seurat", assay = "RNA", overwrite = TRUE)
    
    # "meta.data = FALSE" is to avoid "Error: Missing required datasets 'levels' and 'values'" -> FIXME!!!
    cell_nhood_ann <- paste0(out_miloR_prefix, "_nhoodGroup_annotation.tsv") # by cell
    obj <- LoadH5Seurat(seurat_file, meta.data = FALSE) # still throuws an error, even including misc = FALSE
    df_meta <- read.table(cell_nhood_ann, sep = "\t", header = TRUE)
    if (!is.null(meta_files)) {
        for (meta_file in meta_files) {
            df <- read.table(meta_file, sep = "\t", header = TRUE, check.names = FALSE)[,2,drop=FALSE]
            df_meta <- cbind(df_meta, df)
        }
    }
    obj <- CreateSeuratObject(counts = obj@assays$RNA@counts, meta.data = df_meta)
    obj <- NormalizeData(obj) # this is needed for FC calculation
    
    return(obj)
}


# the object should contain *all* genes 
# the meta_file_list is a vector of file names that contain metadata for the count matrix
# ideally, it should be: nhood_group metadata (cell assignment to nhood_groups), hard cell-type annotation, hard category annotation
annotate_object <- function(adata_prefix, out_miloR_prefix, meta_files = NULL) {

    adata_file <- paste0(adata_prefix, ".h5ad")
    seurat_file <- paste0(adata_prefix, ".h5seurat")
    adata <- read_h5ad(adata_file, to = "InMemoryAnnData")
    
    # FIXED: I only need counts from adata, so use the below workaround (since LoadH5Seurat still throws an error)
    rownames(adata$X) <- adata$obs_names
    colnames(adata$X) <- adata$var_names
    M <- t(adata$X)
    
    # "meta.data = FALSE" is to avoid "Error: Missing required datasets 'levels' and 'values'"
    cell_nhood_ann <- paste0(out_miloR_prefix, "_nhoodGroup_annotation.tsv") # by cell
    df_meta <- read.table(cell_nhood_ann, sep = "\t", header = TRUE)
    if (!is.null(meta_files)) {
        for (meta_file in meta_files) {
            df <- read.table(meta_file, sep = "\t", header = TRUE, check.names = FALSE)[,2,drop=FALSE]
            df_meta <- cbind(df_meta, df)
        }
    }
    
    obj <- CreateSeuratObject(counts = M, meta.data = df_meta)
    obj <- NormalizeData(obj) # this is needed for FC calculation
    
    return(obj)
}


# meta can be obj@meta.data
# use color.by = "cell_type" or color.by = "category"
plot_nhood_group_annotation <- function(meta, da_prefix, out_miloR_prefix, color.by = "cell_type") {

    da_results_nhoodgroup_file <- paste0(out_miloR_prefix, "_withNhoodGroups.tsv") # by nhood
    da.res <- read.table(da_results_nhoodgroup_file, sep = "\t", header = TRUE)
    fields <- grep("NhoodGroup", colnames(meta), value = TRUE)
    for (field in fields) {
        groups <- unique(da.res[[field]])
        groups <- groups[!is.na(groups)]
        avg_logfc <- unlist(lapply(groups, function(x) mean(da.res$logFC[da.res[[field]] == x], na.rm = TRUE)))
        groups <- groups[order(avg_logfc)]
    
        df <- meta
        colnames(df)[colnames(df) == field] <- "NhoodGroup"
        colnames(df)[colnames(df) == color.by] <- "color.by"
        df$NhoodGroup <- factor(df$NhoodGroup, levels = groups)

        # plot cell type distribution in nhood groups
        pdf(paste0(out_miloR_prefix, "_", field, "_barplot_", color.by, ".pdf"), width = 10)
        print(ggplot(data = df[!is.na(df$NhoodGroup),], aes(x = NhoodGroup, fill = color.by)) + geom_bar() + coord_flip()) + ggtitle(field)
        dev.off()
    }
}


find_all_nhood_group_markers <- function(obj, out_file) {
    fields <- grep("NhoodGroup", colnames(obj@meta.data), value = TRUE)
    df_all <- as.data.frame(matrix(NA, nrow = 0, ncol = 6)) # fc, pct.1, pct.2, NhoodGroup, gene, NhoodGroupParams
    for (field in fields) {
        df <- find_nhood_group_markers(obj = obj, meta.field = field)
        df <- cbind(df, NhoodGroupParams = rep(field, nrow(df)))
        rownames(df) <- paste(field, rownames(df), sep = "_")
        df_all = rbind(df_all, df)
    }
    write.table(df_all, file = out_file, sep = "\t", quote = FALSE)
}


# output example (function FoldChange):
#       avg_log2FC pct.1 pct.2
# ISG15   0.4867800 0.645 0.418
# UBE2J2 -0.2808274 0.111 0.161
find_nhood_group_markers <- function(obj, meta.field, min.cells = 10, min.pct = 0.1, min.cells.group = 10, min.logfc = 0.25) {
    Idents(obj) <- meta.field
#    dea <- FindAllMarkers(obj, test.use = test.use) # do not use statistical test here - perhaps on specific nhoods or nhood groups afterwards
    groups <- unique(Idents(obj))
    groups <- groups[!is.na(groups)]
    ncells <- rowSums(obj@assays$RNA@counts > 0)
    features <- rownames(obj)[ncells >= min.cells] # only compute logFC for features detected in >= min.cells cells
    dea <- as.data.frame(matrix(NA, nrow = 0, ncol = 5)) # fc, pct.1, pct.2, cluster, gene
    for (group in groups) {
        idx <- which(Idents(obj) == group)
        if (length(idx) > min.cells.group) { # only compute logFC for groups with size >= min.cells.group
            cells.1 <- colnames(obj)[idx]
            cells.2 <- colnames(obj)[-idx]
            df <- FoldChange(obj, features = features, cells.1, cells.2) # potential issue https://github.com/satijalab/seurat/issues/6701 -> OK
            idx2 <- which(df$pct.1 >= min.pct & df$pct.2 >= min.pct & abs(df$avg_log2FC) >= min.logfc) 
            if (length(idx2) > 0) {
                df <- df[idx2,]
                gene = rownames(df)
                NhoodGroup = rep(group, nrow(df))
                rownames(df) <- paste(group, gene, sep = "_")
                dea <- rbind(dea, cbind(df, gene, NhoodGroup))
            }
        }
    }
    return(dea)
}


# OLD FUNCTION: superseded by the above functions
group_neighbourhoods <- function(sce_file, da_results_file, out_prefix, cat_file, ct_file) {

    sce <- readRDS(sce_file)
    
    # create cellxnhoodgroup matrix
    # sce@nhoods: rows = cells, columns = nhoods
    # da.res$X: nhood ID
    # da.res$NhoodGroup: nhood group
    groups <- unique(da.res$NhoodGroup)
    M1 <- sapply(groups, function(x) as.numeric(da.res$NhoodGroup == x))
    M1 <- Matrix(M1) # nhoodxnhoodgroup matrix
    rownames(M1) <- da.res$X
    colnames(M1) <- groups
    M <- sce@nhoods %*% M1 # cellxnhoodgroup matrix containing the count of nhood
    
    # assign each cell to a group by selecting the one with highest number of nhood containing the cell
    # delete the assignment if the top groups contain less than MIN_PERC fraction of nhoods that contain the cell
    assignment <- apply(M, 1, which.max)
    max.val <- apply(M/rowSums(M), 1, max)
    idx <- which(max.val < MIN_PERC)
    assignment[idx] <- rep(NaN, length(idx))
    
    avg_logfc <- unlist(lapply(groups, function(x) mean(da.res$logFC[da.res$NhoodGroup == x])))
    groups <- groups[order(avg_logfc)]
    
    # create metadata
    ct <- read.table(ct_file, sep = "\t", header = TRUE)[,2]
    cat <- read.table(cat_file, sep = "\t", header = TRUE)[,2]
    meta <- data.frame(condition = sce$condition, NhoodGroup = assignment, cell.type = ct, category = cat)
    meta$NhoodGroup <- factor(meta$NhoodGroup, levels = groups)

    # plot cell type distribution in nhood groups
    pdf(paste0(out_prefix, "_NhoodGroup_barplot_cell_type.pdf"), width = 10)
    print(ggplot(data = meta[!is.na(meta$NhoodGroup),], aes(x = NhoodGroup, fill = cell.type)) + geom_bar() + coord_flip())
    dev.off()

    # plot category distribution in nhood groups
    pdf(paste0(out_prefix, "_NhoodGroup_barplot_category.pdf"))
    print(ggplot(data = meta[!is.na(meta$NhoodGroup),], aes(x = NhoodGroup, fill = category)) + geom_bar() + coord_flip())
    dev.off()
    
}



