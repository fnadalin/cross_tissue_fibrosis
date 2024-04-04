
# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

library("anndataR")
library("SingleCellExperiment")
library("miloR")
library("dplyr")

set.seed(123)

# FOR TESTING: use the atlas sample dataset with "disease" assigned to the last 500 cells: tests/atlas_liver_sample_simConditions.h5ad
# (in the query sample dataset, no DA nhood was found)
# N.B.: make sure that cells in the same sample are all assigned to the same condition!!!
differential_abundance_milo <- function(adata_file, latent_id = 'X_scvi', out_prefix) {
    
    # convert adata to SingleCellExperiment
    adata <- read_h5ad(adata_file, to = "InMemoryAnnData")
    sce <- adata$to_SingleCellExperiment()
    
    # prepare the SingleCellExperiment object for miloDE
    colnames(sce) <- adata$obs$barcode
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
# TODO: maybe reformat labels_hard as a 0-1 matrix, to be consistent with the soft label format
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


