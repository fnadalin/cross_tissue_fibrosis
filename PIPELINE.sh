
source ~/.bashrc # for conda init

# PARAMS="params/param_Set1.tsv"

if [[ $# -lt 1 ]]
then
    echo "Usage: $0 <param_file>"
    exit
fi

PARAMS=$1

DIR=$( cd $(dirname $0) ; pwd )

create embeddings
conda activate scvi-env
python ${DIR}/1_latent_embedding.py ${PARAMS}
conda deactivate

# differential abundance with MiloR
# N.B. sample_id is hardcoded here...
module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf # NEW!!!
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER
Rscript ${DIR}/2_differential_abundance.R ${PARAMS}

# differential expression with MiloDE
R_LIBS_USER=""
export R_LIBS_USER
singularity exec /nfs/research/marioni/andrian/containers/miloDE_cms.simg Rscript ${DIR}/2_differential_expression.R ${PARAMS}

# differential analysis VS cell type annotation
# conda activate sklearn-env
# python ${DIR}/3_multinomial_regression.py ${PARAMS}
# conda deactivate

exit

