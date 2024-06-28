
DIR <- "src/cross_tissue_fibrosis/"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <in_dir> <out_dir> <adata>\n")
    q()
}

in_dir <- args[1]
out_dir <- args[2]
cells <- args[3]

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
out_prefix <- file.path(out_dir, "Intestine_Xavier")
load_Xavier(in_dir, out_prefix, cells)

q()

