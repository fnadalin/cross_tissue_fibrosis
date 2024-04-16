
DIR <- "src/cross_tissue_fibrosis/"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <cell_obj> <nuclei_obj> <out_prefix>\n")
    q()
}

cell_obj <- args[1]
nuclei_obj <- args[2]
out_dir <- args[3]

### export functions

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
source(file.path(DIR, "datasets.R"))
setwd(WORKING_DIR)

### execute

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cell_out_prefix <- file.path(out_dir, "Heart_Lavine_cell")
nuclei_out_prefix <- file.path(out_dir, "Heart_Lavine_nuclei")
seurat_to_adata_Lavine_cells(cell_obj, cell_out_prefix)
seurat_to_adata_Lavine_nuclei(nuclei_obj, nuclei_out_prefix)

q()

