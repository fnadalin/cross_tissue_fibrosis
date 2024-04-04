
PARAMS="params/param_Set1.tsv"

# create embeddings
conda activate scvi-env
python 1_latent_embedding.py ${PARAMS}
conda deactivate

# differential abundance with MiloR
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER
Rscript 2_differential_abundance.R ${PARAMS}

# differential expression with MiloDE
R_LIBS_USER=""
export R_LIBS_USER
singularity exec /nfs/research/marioni/andrian/containers/miloDE_cms.simg Rscript 2_differential_expression.R ${PARAMS}

# differential analysis VS cell type annotation
conda activate sklearn-env
python 3_multinomial_regression.py ${PARAMS}
conda deactivate

exit

