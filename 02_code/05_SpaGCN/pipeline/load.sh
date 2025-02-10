##!/usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --matrixPath)
            matrixPath="$2"
            shift # past argument
            shift # past value
            ;;
        --tissuePositionPath)
            tissuePositionPath="$2"
            shift # past argument
            shift # past value
            ;;
        --imgPath)
            imgPath="$2"
            shift # past argument
            shift # past value
            ;;
        --outPath)
            outPath="$2"
            shift # past argument
            shift # past value
            ;;
        --mainscriptsFolder)
            mainscriptsFolder="$2"
            shift # past argument
            shift # past value
            ;;
        --parentJobName)
            parentJobName="$2"
            shift # past argument
            shift # past value
            ;;
        --alphaArray)
            alphaArray="$2"
            shift # past argument
            shift # past value
            ;;
        --betaArray)
            betaArray="$2"
            shift # past argument
            shift # past value
            ;;
        --pArray)
            pArray="$2"
            shift # past argument
            shift # past value
            ;;
        --nClustersArray)
            nClustersArray="$2"
            shift # past argument
            shift # past value
            ;;
        --toLoad)
            toLoad="$2"
            shift # past argument
            shift # past value
            ;;
        --toAdj)
            toAdj="$2"
            shift # past argument
            shift # past value
            ;;
        --toDetectDomain)
            toDetectDomain="$2"
            shift # past argument
            shift # past value
            ;;
        --toDetectSVG)
            toDetectSVG="$2"
            shift # past argument
            shift # past value
            ;;
        --toCallBack)
            toCallBack="$2"
            shift # past argument
            shift # past value
            ;;
        --callBackPath)
            callBackPath="$2"
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

echo "run parameters: \
matrixPath=${matrixPath} \
tissuePositionPath=${tissuePositionPath} \
imgPath=${imgPath} \
outPath=${outPath} \
mainscriptsFolder=${mainscriptsFolder} \
parentJobName=${parentJobName} \
alphaArray=${alphaArray} \
betaArray=${betaArray} \
pArray=${pArray} \
nClustersArray=${nClustersArray} \
toLoad=${toLoad} \
toAdj=${toAdj} \
toDetectDomain=${toDetectDomain} \
toDetectSVG=${toDetectSVG} \
toCallBack=${toCallBack} \
callBackPath=${callBackPath}
"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

export CV_IO_MAX_IMAGE_PIXELS=1099511627776
source ${HOME}/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate SpaGCN


rawAdataPath=${outPath}/raw_adata.h5ad
if [ "$toLoad" == "YES" ]; then
    python /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/src/step1_load_adata.py \
        --matrixPath=${matrixPath} \
        --tissuePositionPath=${tissuePositionPath} \
        --outPath=${rawAdataPath}
fi

if [ "$toAdj" == "NO" ]; then
    mmForInData=50
    qName="e80short"
    wTime="2:50"
else
    mmForInData=400
    qName="e80medium"
    wTime="23:50"
fi


IFS='-' read -a alphaArray <<< "$alphaArray"
IFS='-' read -a betaArray <<< "$betaArray"

for alpha in "${alphaArray[@]}"; do
    for beta in "${betaArray[@]}"; do
        tempFolder=${outPath}/alpha_${alpha}_beta_${beta}
        if [ ! -d $tempFolder ]; then
            mkdir -p $tempFolder
        fi
        
        JOBFOLDER=${tempFolder}
        JOBNAME=${parentJobName}_load_alpha_${alpha}_beta_${beta}
        if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
            rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
        fi

        bsub \
            -J ${JOBNAME} \
            -o ${JOBFOLDER}/${JOBNAME}.o.txt \
            -e ${JOBFOLDER}/${JOBNAME}.e.txt \
            -cwd ${JOBFOLDER} \
            -q ${qName} \
            -W ${wTime} \
            -n 1 \
            -M ${mmForInData} \
            -R rusage[mem=${mmForInData}] \
            /bin/bash -c " \
                /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/pipeline/adj.sh \
                --adataPath ${rawAdataPath} \
                --imgPath ${imgPath} \
                --outPath ${JOBFOLDER} \
                --mainscriptsFolder ${mainscriptsFolder} \
                --parentJobName ${parentJobName} \
                --alpha ${alpha} \
                --beta ${beta} \
                --pArray ${pArray} \
                --nClustersArray ${nClustersArray} \
                --toAdj ${toAdj} \
                --toDetectDomain ${toDetectDomain} \
                --toDetectSVG ${toDetectSVG} \
                --toCallBack ${toCallBack} \
                --callBackPath ${callBackPath} "
    done
done
