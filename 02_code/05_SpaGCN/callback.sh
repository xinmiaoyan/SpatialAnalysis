##!/usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --adataPath)
            adataPath="$2"
            shift # past argument
            shift # past value
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "run parameters: adataPath=${adataPath}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate SpaGCN
export CV_IO_MAX_IMAGE_PIXELS=1099511627776


/rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/run_SpaGCN/src/step4_module_score_plot.py \
    --adataPath=${adataPath} \
    --markerTablePath="/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/ImmuneCellMarkers.tsv"
