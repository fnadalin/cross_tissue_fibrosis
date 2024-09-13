library("Seurat")
library("SeuratDisk")
library("Matrix")

Kaminski_genes <- "GSE136831_AllCells.GeneIDs.txt"
Kaminski_meta <- "GSE136831_AllCells.Samples.CellType.MetadataTable.txt"
Kaminski_barcodes <- "GSE136831_AllCells.cellBarcodes.txt"
Kaminski_matrix <- "GSE136831_RawCounts_Sparse.mtx.gz"

Xavier_cond_CO <- c("CO_EPI", "CO_IMM", "CO_STR")
Xavier_cond_TI <- c("TI_EPI", "TI_IMM", "TI_STR")
Xavier_meta <- "scp_metadata_combined.v2.txt"

Helmsley_normal <- c()
Helmsley_disease <- c() 


# FIXME: check the output objects
# FIXME: check the metadata / other info in Lavine because there is inconsistency in # of features!!!

seurat_to_adata_Lavine_cells <- function(obj, out_prefix) {
    load(obj)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
#    HDCM@active.assay = "RNA"
#    SaveH5Seurat(HDCM, filename = out_seurat, overwrite = TRUE)
    #### NEW ####
    meta <- HDCM@meta.data
    meta <- cbind(meta, assay = rep("cell", nrow(meta)))
    object <- CreateSeuratObject(counts = HDCM@assays$RNA@counts, project = "Lavine", meta.data = meta)
    SaveH5Seurat(object, filename = out_seurat, overwrite = TRUE)
    #############
    Convert(out_seurat, dest = "h5ad", overwrite = TRUE)
}

seurat_to_adata_Lavine_nuclei <- function(obj, out_prefix) {
    load(obj)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
#    nuclei@active.assay = "RNA"
#    SaveH5Seurat(nuclei, filename = out_seurat, overwrite = TRUE)
    #### NEW ####
    meta <- nuclei@meta.data
    meta <- cbind(meta, assay = rep("nuclei", nrow(meta)))
    object <- CreateSeuratObject(counts = nuclei@assays$RNA@counts, project = "Lavine", meta.data = meta)
    SaveH5Seurat(object, filename = out_seurat, overwrite = TRUE)
    #############
    Convert(out_seurat, dest = "h5ad", overwrite = TRUE)
}

seurat_to_adata_Helmsley <- function(obj, out_prefix, condition) {
    obj <- readRDS(obj)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
#    HDCM@active.assay = "RNA"
#    SaveH5Seurat(HDCM, filename = out_seurat, overwrite = TRUE)
    #### NEW ####
    cells <- colnames(obj)[!is.na(obj@meta.data$cell_type_final)]
    obj <- subset(obj, cells = cells)
    meta <- obj@meta.data
    meta <- cbind(meta, condition = rep(condition, nrow(meta)))
    object <- CreateSeuratObject(counts = obj@assays$RNA@counts, project = "Helmsley", meta.data = meta)
    SaveH5Seurat(object, filename = out_seurat, overwrite = TRUE)
    #############
    Convert(out_seurat, dest = "h5ad", overwrite = TRUE)
}

load_Kaminski <- function(in_dir, out_prefix) {
    matrix_file <- file.path(in_dir, Kaminski_matrix)
    genes_file <- file.path(in_dir, Kaminski_genes)
    barcode_file <- file.path(in_dir, Kaminski_barcodes)
    metadata_file <- file.path(in_dir, Kaminski_meta)
    
    M <- readMM(gzfile(matrix_file))
    genes = read.table(genes_file, sep = "\t", header = TRUE) 
    barcodes = read.table(barcode_file, sep = "\t", header = FALSE)[,1]
    meta = read.table(metadata_file, sep = "\t", header = TRUE)
    colnames(M) <- rownames(meta) <- barcodes
    rownames(M) <- genes[['HGNC_EnsemblAlt_GeneID']]
    
    obj <- CreateSeuratObject(counts = M, project = "Kaminski", meta.data = meta)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
    SaveH5Seurat(obj, filename = out_seurat, overwrite = TRUE)
    Convert(out_seurat, dest = "h5ad", overwrite = TRUE)
}

load_Xavier <- function(in_dir, out_prefix, cells_file) {
    load_Xavier_cond(in_dir, paste0(out_prefix, "_CO"), cells_file, conds = Xavier_cond_CO)
    load_Xavier_cond(in_dir, paste0(out_prefix, "_TI"), cells_file, conds = Xavier_cond_TI)
}

load_Xavier_cond <- function(in_dir, out_prefix, cells_file, conds) {
    feature_file <- file.path(in_dir, paste0(conds[1], ".scp.features.tsv"))
    features <- read.table(feature_file, sep = "\t")
    matrix <- Matrix(nrow = nrow(features), ncol = 0)
    rownames(matrix) <- features[,2] # keep the gene symbol
    condition <- NULL
    for (cond in conds) {
        barcode_file <- file.path(in_dir, paste0(cond, ".scp.barcodes.tsv"))
        matrix_file <- file.path(in_dir, paste0(cond, ".scp.raw.mtx"))
        matrix_tmp <- readMM(matrix_file)
        colnames(matrix_tmp) <- read.table(barcode_file, sep = "\t")[,1]
        matrix <- cbind2(matrix, matrix_tmp)
        condition <- c(condition, rep(cond, ncol(matrix_tmp)))
    }
    
    meta_file <- file.path(in_dir, Xavier_meta)
    meta <- read.table(meta_file, sep = "\t", header = TRUE)
    meta <- meta[2:nrow(meta),]
    cells <- read.table(cells_file, sep = "\t", header = TRUE)[,2]
    meta <- meta[meta[,"NAME"] %in% cells,]
    cells <- colnames(matrix)[colnames(matrix) %in% meta[,"NAME"]]
    idx <- colnames(matrix) %in% cells
    matrix <- matrix[,idx]
    condition <- condition[idx]
    meta <- meta[meta[,"NAME"] %in% cells,]
    meta <- meta[match(cells, meta[,"NAME"]),]
    meta <- cbind(meta, condition = condition)
    
    obj <- CreateSeuratObject(counts = matrix, project = "Xavier", meta.data = meta)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
    SaveH5Seurat(obj, filename = out_seurat, overwrite = TRUE)
    Convert(out_seurat, dest = "h5ad", overwrite = TRUE)
}


