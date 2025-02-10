#BSUB -J submit
#BSUB -q highmem
#BSUB -W 220:00
#BSUB -n 8
#BSUB -M 1000
#BSUB -R rusage[mem=1000]
#BSUB -B
#BSUB -N
#BSUB -u xyan4@mdanderson.org
#BSUB -o /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/01_spaceranger_mkfastq/submit.o.txt
#BSUB -e /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/01_spaceranger_mkfastq/submit.e.txt
#BSUB -cwd //rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/01_spaceranger_mkfastq/
rm -rf /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/01_spaceranger_mkfastq/submit.o.txt
rm -rf /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/01_spaceranger_mkfastq/submit.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____

PROJECT_FOLDER=/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng
DATA_FOLDER=${PROJECT_FOLDER}/01_data
RESULT_FOLDER=${PROJECT_FOLDER}/03_result
CODE_FOLDER=${PROJECT_FOLDER}/02_code
PIPELINE_NAME=01_spaceranger_mkfastq
PROJECT_NAME=$(basename ${PROJECT_FOLDER})

OutDir=$RESULT_FOLDER/$PIPELINE_PATH_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi


cd $OutDir
module load bcl2fastq spaceranger

spaceranger mkfastq --id=BHK2N3DMXY --run=/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/230509_A00422_0881_BH75LKDRX3 --csv=/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/230509_A00422_0881_BH75LKDRX3/SampleSheet10X.csv





