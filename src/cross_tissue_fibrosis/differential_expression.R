# for converting from/to anndata/singlecellexperiment objects, see https://github.com/scverse/anndataR

# singularity exec /nfs/research/marioni/andrian/containers/miloDE_cms.simg Rscript
library("anndataR", lib.loc = "/hps/software/users/marioni/francesca/R_libs")
library("SingleCellExperiment")
library("miloDE")
library("miloR")

set.seed(123)

# FOR TESTING: use the atlas sample dataset with "disease" assigned to the last 500 cells: tests/atlas_liver_sample_simConditions.h5ad
# (in the query sample dataset, no DA nhood was found) 
differential_expression_miloDE <- function(adata_file, latent_id = 'X_scvi', out_prefix = "tests/miloDE/out") {
    
    # convert adata to SingleCellExperiment
    adata <- read_h5ad(adata_file, to = "InMemoryAnnData")
    sce <- adata$to_SingleCellExperiment()
    
    # prepare the SingleCellExperiment object for miloDE
    colnames(sce) <- adata$obs_names # NEW: assign cell names
    names(assays(sce)) <- "counts"
    reducedDim(sce) <- adata$obsm[[latent_id]]
    reducedDimNames(sce) <- latent_id

    sce <- assign_neighbourhoods(sce, k = 20, order = 2, filtering = TRUE, reducedDim_name = latent_id)
    # the above gives the following error (on the git version):
    # Filtering redundant neighbourhoods.
    # Error in if (isEmpty(nhoodIndex(x))) { : the condition has length > 1 ---> even in the test dataset!!!
    de_stat <- de_test_neighbourhoods(sce, sample_id = "sample_id", design = ~condition, covariates = c("condition"))

    out_de <- paste0(out_prefix, ".tsv")
    write.table(de_stat, file = out_de, quote = FALSE, sep = "\t", row.names = FALSE)
    
    out_sce <- paste0(out_prefix, ".Rds")
    saveRDS(sce, file = out_sce)

    # the below code doesn't work...
#    miloDE_out <- adata$uns
#    miloDE_out["nhoods"] <- sce@nhoods
#    miloDE_out["nhoodIndex"] <- NULL
#    append(miloDE_out["nhoodIndex"], sce@nhoodIndex)
#    miloDE_out["nhoodExpression"] <- sce@nhoodExpression
#    miloDE_out["nhoodAdjacency"] <- sce@nhoodAdjacency
#    adata$uns["miloDE_out"] <- milo_out
    
#    out_adata <- paste0(out_prefix, ".h5ad")
#    write_h5ad(adata, path = out_adata)
    
}


# adapted from annotateNhoods from package:miloR
# compute a nhood label probability by majority voting
# TODO: maybe reformat labels_hard as a 0-1 matrix, to be consistent with the soft label format
annotate_neighbourhoods_soft_miloDE <- function(sce_file, de_results_file, soft_ann_file) {

    x <- readRDS(sce_file)
    de.res <- read.table(de_results_file, sep = "\t", header = TRUE)
    soft_ann <- read.table(soft_ann_file, sep = "\t", header = TRUE, row.names = 1)

    if (!is(x, "Milo")) {
        stop("Unrecognised input type - must be of class Milo")
    }
    if (ncol(nhoods(x)) != length(unique(de.res$Nhood))) {
        stop("the number of rows in de.res does not match the number of neighbourhoods in nhoods(x). Are you sure de.res is the output of de_test_neighbourhoods(x)?")
    }
    if (nrow(nhoods(x)) != nrow(soft_ann)) {
        stop("the number of rows in soft_ann does not match the number of cells in nhoods(x). Are you sure Milo was run on the same object?")
    }

    nhood_counts <- vapply(seq_len(ncol(nhoods(x))), FUN = function(n) colSums(soft_ann[which(nhoods(x)[,n] == 1),]), FUN.VALUE = numeric(ncol(soft_ann)))
    nhood_counts <- nhood_counts/colSums(nhood_counts)
    nhood_counts <- t(nhood_counts)
    rownames(nhood_counts) <- seq_len(ncol(nhoods(x)))
    colnames(nhood_counts) <- colnames(soft_ann)
    de.res <- cbind(de.res, nhood_counts[de.res$Nhood,])
    
    return(de.res)
}


write_annotate_neighbourhoods_soft_miloDE <- function(de.res, de_results_ann_file) {

    df <- data.frame(rownames(de.res))
    de.res <- cbind(df, de.res)
    colnames(de.res)[1] <- ""
    write.table(de.res, file = de_results_ann_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


# evaluate whether another method is needed...
differential_expression_singleCellHaystack <- function(adata_file, latent_id = 'X_scVI', design = "~ disease", out_file) {}


