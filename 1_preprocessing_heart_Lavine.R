
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <cell_obj> <nuclei_obj> <out_prefix>\n")
    q()
}

cell_obj <- args[1]
nuclei_obj <- args[2]
out_prefix <- args[3]

cell_out_prefix <- file.path(out_prefix, "Heart_Lavine_cell")
nuclei_out_prefix <- file.path(out_prefix, "Heart_Lavine_nuclei")
seurat_to_adata_Lavine_cells(cell_obj, cell_out_prefix)
seurat_to_adata_Lavine_nuclei(nuclei_obj, nuclei_out_prefix)

q()

