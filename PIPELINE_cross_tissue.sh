
if [[ $# -lt 2 ]]
then
    echo "Usage: $0 <param1_file> <param2_file>"
    exit
fi

PARAMS1=$1
PARAMS2=$2

DIR=$( cd $(dirname $0) ; pwd )

# cross-tissue DA nhood group comparison
# Rscript ${DIR}/nhg_comparison_across_tissues.R ${PARAMS1} ${PARAMS2} # OLD -> not run!
Rscript ${DIR}/1_nhg_comparison_across_tissues.R ${PARAMS1} 
Rscript ${DIR}/1_nhg_comparison_across_tissues.R ${PARAMS2} 
Rscript ${DIR}/2_nhg_comparison_across_tissues.R ${PARAMS1} ${PARAMS2}

exit

