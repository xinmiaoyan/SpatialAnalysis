#BSUB -J submit
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/submit.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/submit.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/submit.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/submit.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____

PROJECT_FOLDER=/rsrch3/scratch/genomic_med/ychu2/projects/project26
DATA_FOLDER=${PROJECT_FOLDER}/data
RESULT_FOLDER=${PROJECT_FOLDER}/result
CODE_FOLDER=${PROJECT_FOLDER}/code
PIPELINE_FOLDER=${CODE_FOLDER}/pipeline
SRC_FOLDER=${CODE_FOLDER}/src
KNOWLEDGE_FOLDER=${PROJECT_FOLDER}/knowledge
PIPELINE_NAME=PM3131__3_SpaGCN
PIPELINE_PATH_NAME=PM3131/3_SpaGCN
PROJECT_NAME=$(basename ${PROJECT_FOLDER})

OutDir=$RESULT_FOLDER/$PIPELINE_PATH_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi


image_Array=(
    /rsrch3/scratch/genomic_med/ychu2/projects/project26/data/PM3131/img/HandE/HKU-01_15BX24742.tif
    /rsrch3/scratch/genomic_med/ychu2/projects/project26/data/PM3131/img/HandE/HKU-02_16BX12829.tif
    /rsrch3/scratch/genomic_med/ychu2/projects/project26/data/PM3131/img/HandE/HKU-03_16BX16714.tif
    /rsrch3/scratch/genomic_med/ychu2/projects/project26/data/PM3131/img/HandE/HKU-04_16BX28164.tif )

sample_folder_Array=( HKU01 HKU02 HKU03 HKU04 )

for i in "${!sample_folder_Array[@]}"; do 
    sample_folder=${sample_folder_Array[$i]}
    matrixPath=/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/2_spaceranger_count/${sample_folder}/outs/filtered_feature_bc_matrix.h5
    tissuePositionPath=/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/2_spaceranger_count/${sample_folder}/outs/spatial/tissue_positions.csv
    imgPath=${image_Array[$i]}
    outPath=${OutDir}/${sample_folder}
    if [ ! -d $outPath ]; then
        mkdir -p $outPath
    fi
    mainscriptsFolder=$(pwd)
    parentJobName=SpaGCN
    alphaArray="0.5-1-1.5"
    betaArray="30-49-70"
    pArray="0.4-0.5-0.6"
    nClustersArray="6-7-8-9-10-11-12-13-14-15"
    # alphaArray="1"
    # betaArray="49"
    # pArray="0.5"
    # nClustersArray="3"
    toLoad="NO"
    toAdj="NO"
    toDetectDomain="NO"
    toDetectSVG="NO"
    toCallBack="YES"
    callBackPath=/rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/callback.sh

    /rsrch3/scratch/genomic_med/ychu2/projects/project26/code/pipeline/private/PM3131/3_SpaGCN/run_SpaGCN/pipeline/run.sh \
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
        --callBackPath ${callBackPath}

done
