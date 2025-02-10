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
        --alpha)
            alpha="$2"
            shift # past argument
            shift # past value
            ;;
        --beta)
            beta="$2"
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

echo "run parameters:
adataPath=${adataPath}
imgPath=${imgPath}
outPath=${outPath}
mainscriptsFolder=${mainscriptsFolder}
parentJobName=${parentJobName}
alpha=${alpha}
beta=${beta}
toAdj=${toAdj}
toDetectDomain=${toDetectDomain}
toDetectSVG=${toDetectSVG}
toCallBack=${toCallBack}
callBackPath=${callBackPath}
"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi


echo "$(which python)"
echo $(which python)

export CV_IO_MAX_IMAGE_PIXELS=1099511627776
source ${HOME}/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate SpaGCN

adjPath=${outPath}/adj.tsv
if [ "$toAdj" == "YES" ]; then
    echo "$(which python)"
    echo $(which python)
    echo "$(python --version)"
    echo $(python --version)
    python /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/src/step2_cal_adj_matrix.py \
        --adataPath=${adataPath} \
        --imgPath=${imgPath} \
        --outPath=${adjPath} \
        --alpha=${alpha} \
        --beta=${beta}
fi

if [ "$toDetectDomain" == "NO" ]; then
    mmForInData=50
    qName="e80short"
    wTime="2:50"
else
    mmForInData=400
    qName="e80long"
    wTime="47:50"
fi


IFS='-' read -a pArray <<< "$pArray"
IFS='-' read -a nClustersArray <<< "$nClustersArray"
for p in "${pArray[@]}"; do
    for nClusters in "${nClustersArray[@]}"; do
        tempFolder=${outPath}/p_${p}_nClusters_${nClusters}
        if [ ! -d $tempFolder ]; then
            mkdir -p $tempFolder
        fi
        
        JOBFOLDER=${tempFolder}
        JOBNAME=${parentJobName}_load_p_${p}_nClusters_${nClusters}
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
                /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/pipeline/domain.sh \
                --adataPath ${adataPath} \
                --outPath ${JOBFOLDER} \
                --adjPath ${adjPath} \
                --mainscriptsFolder ${mainscriptsFolder} \
                --parentJobName ${parentJobName} \
                --p ${p} \
                --nClusters ${nClusters} \
                --toDetectDomain ${toDetectDomain} \
                --toDetectSVG ${toDetectSVG} \
                --toCallBack ${toCallBack} \
                --callBackPath ${callBackPath} "
    done
done

