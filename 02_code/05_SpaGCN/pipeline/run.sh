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

echo "run parameters:
matrixPath=${matrixPath}
tissuePositionPath=${tissuePositionPath}
imgPath=${imgPath}
outPath=${outPath}
mainscriptsFolder=${mainscriptsFolder}
parentJobName=${parentJobName}
alphaArray=${alphaArray}
betaArray=${betaArray}
pArray=${pArray}
nClustersArray=${nClustersArray}
toLoad=${toLoad}
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

# source ${HOME}/miniconda3/etc/profile.d/conda.sh
# conda activate SpaGCN

if [ "$toLoad" == "NO" ]; then
    mmForInData=50
    qName="e80short"
    wTime="2:50"
else
    mmForInData=400
    qName="e80medium"
    wTime="23:50"
fi

if [ ! -d $outPath ]; then
    mkdir -p $outPath
fi

JOBFOLDER=${outPath}
JOBNAME=${parentJobName}_load
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
        /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/pipeline/load.sh \
            --matrixPath ${matrixPath} \
            --tissuePositionPath ${tissuePositionPath} \
            --imgPath ${imgPath} \
            --outPath ${outPath} \
            --mainscriptsFolder ${mainscriptsFolder} \
            --parentJobName ${parentJobName} \
            --alphaArray ${alphaArray} \
            --betaArray ${betaArray} \
            --pArray ${pArray} \
            --nClustersArray ${nClustersArray} \
            --toLoad ${toLoad} \
            --toAdj ${toAdj} \
            --toDetectDomain ${toDetectDomain} \
            --toDetectSVG ${toDetectSVG} \
            --toCallBack ${toCallBack} \
            --callBackPath ${callBackPath} "

echo  " \
        /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/pipeline/load.sh \
            --matrixPath ${matrixPath} \
            --tissuePositionPath ${tissuePositionPath} \
            --imgPath ${imgPath} \
            --outPath ${outPath} \
            --mainscriptsFolder ${mainscriptsFolder} \
            --parentJobName ${parentJobName} \
            --alphaArray ${alphaArray} \
            --betaArray ${betaArray} \
            --pArray ${pArray} \
            --nClustersArray ${nClustersArray} \
            --toLoad ${toLoad} \
            --toAdj ${toAdj} \
            --toDetectDomain ${toDetectDomain} \
            --toDetectSVG ${toDetectSVG} \
            --toCallBack ${toCallBack} \
            --callBackPath ${callBackPath} "
