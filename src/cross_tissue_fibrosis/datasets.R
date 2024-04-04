library("Seurat")
library("SeuratDisk")
library("Matrix")

Kaminski_genes <- "GSE136831_AllCells.GeneIDs.txt"
Kaminski_meta <- "GSE136831_AllCells.Samples.CellType.MetadataTable.txt"
Kaminski_barcodes <- "GSE136831_AllCells.cellBarcodes.txt"
Kaminski_matrix <- "GSE136831_RawCounts_Sparse.mtx.gz"

# FIXME: check the output objects

seurat_to_adata_Lavine_cells <- function(obj, out_prefix) {
    load(obj)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
    SaveH5Seurat(HDCM, filename = out_seurat)
    Convert(out_seurat, dest = "h5ad")
}

seurat_to_adata_Lavine_nuclei <- function(obj, out_prefix) {
    load(obj)
    out_seurat <- paste0(out_prefix, ".h5Seurat")
    SaveH5Seurat(nuclei, filename = out_seurat)
    Convert(out_seurat, dest = "h5ad")
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
    SaveH5Seurat(obj, filename = out_seurat)
    Convert(out_seurat, dest = "h5ad")
}

