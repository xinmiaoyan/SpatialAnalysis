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
        --adjPath)
            adjPath="$2"
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
        --p)
            p="$2"
            shift # past argument
            shift # past value
            ;;
        --nClusters)
            nClusters="$2"
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
adjPath=${adjPath}
outPath=${outPath}
mainscriptsFolder=${mainscriptsFolder}
parentJobName=${parentJobName}
p=${p}
nClusters=${nClusters}
toDetectDomain=${toDetectDomain}
toDetectSVG=${toDetectSVG}
toCallBack=${toCallBack}
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

newAdataPath=${outPath}/domain.h5ad

if [ "$toDetectDomain" == "YES" ]; then
    python /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/src/step3_detect_domain.py \
        --adataPath=${adataPath} \
        --adjPath=${adjPath} \
        --outPath=${newAdataPath}\
        --p=${p} \
        --nClusters=${nClusters}
    python /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/src/step4_plot.py \
        --adataPath=${newAdataPath} \
        --outPath=${outPath}
fi

if [ "$toDetectSVG" == "YES" ]; then
    mmForInData=400
    qName="e80medium"
    wTime="23:50"

    targetArray=($(seq 0 $(( ${nClusters} - 1 ))))

    for target in "${targetArray[@]}"; do
        tempFolder=${outPath}
        if [ ! -d $tempFolder ]; then
            mkdir -p $tempFolder
        fi
        
        JOBFOLDER=${tempFolder}
        JOBNAME=${parentJobName}_SVG_target_${target}
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
                /rsrch3/home/genomic_med/ychu2/configs/public/pipeline/ST/visium/run_SpaGCN/pipeline/SVG.sh \
                    --rawAdataPath ${adataPath} \
                    --adataPath ${newAdataPath} \
                    --outDirPath ${JOBFOLDER} \
                    --target ${target} \
                    --toDetectSVG ${toDetectSVG} "

    done

fi

if [ "$toCallBack" == "YES" ]; then
    mmForInData=10
    qName="e40medium"
    wTime="23:20"

    JOBFOLDER=$(dirname ${newAdataPath})
    JOBNAME=${parentJobName}_CallBack
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
        /bin/bash -c "${callBackPath} --adataPath ${newAdataPath}"
fi
