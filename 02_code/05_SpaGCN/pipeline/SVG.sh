##!/usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --rawAdataPath)
            rawAdataPath="$2"
            shift # past argument
            shift # past value
            ;;
        --adataPath)
            adataPath="$2"
            shift # past argument
            shift # past value
            ;;
        --outDirPath)
            outDirPath="$2"
            shift # past argument
            shift # past value
            ;;
        --target)
            target="$2"
            shift # past argument
            shift # past value
            ;;
        --toDetectSVG)
            toDetectSVG="$2"
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

echo "run parameters:
rawAdataPath=${rawAdataPath}
adataPath=${adataPath}
outDirPath=${outDirPath}
target=${target}
toDetectSVG=${toDetectSVG}
"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

source ${HOME}/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate SpaGCN

export CV_IO_MAX_IMAGE_PIXELS=1099511627776


if [ "$toDetectSVG" == "YES" ]; then
    python /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/src/step4_SVG.py \
        --rawAdataPath=${rawAdataPath} \
        --adataPath=${adataPath} \
        --outDirPath=${outDirPath} \
        --target=${target}
fi
