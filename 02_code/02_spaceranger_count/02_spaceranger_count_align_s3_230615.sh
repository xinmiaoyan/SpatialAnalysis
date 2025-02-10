#BSUB -J submit
#BSUB -q short
#BSUB -W 3:00
#BSUB -n 1
#BSUB -M 50
#BSUB -R rusage[mem=50]
#BSUB -B
#BSUB -N
#BSUB -u xyan4@mdanderson.org
#BSUB -o /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/02_spaceranger_count/submit.o.txt
#BSUB -e /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/02_spaceranger_count/submit.e.txt
#BSUB -cwd /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/02_spaceranger_count/
rm -rf /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/02_spaceranger_count/submit.o.txt
rm -rf /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/02_spaceranger_count/submit.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____

PROJECT_FOLDER=/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng
DATA_FOLDER=${PROJECT_FOLDER}/01_data
RESULT_FOLDER=${PROJECT_FOLDER}/03_result
CODE_FOLDER=${PROJECT_FOLDER}/02_code
PIPELINE_NAME=02_spaceranger_count_S13_S15_aglin
PROJECT_NAME=$(basename ${PIPELINE_NAME})

OutDir=$RESULT_FOLDER/$PROJECT_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi
cd $OutDir

image_Array=(
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/S13-014988_A2.tif
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/S15-061343_A1.tif
    )

cyt_image_Array=(
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/CAVG10036_2023-03-01_17-53-28_chemorun1_V42Y23-325_D1_s13014988.tif
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/CAVG10036_2023-03-01_17-53-28_chemorun1_V42Y23-325_A1_s15061343.tif
    )

slide_Array=(
    V42Y23-325
    V42Y23-325
     )

area_Array=(
    D1
    A1 )

slidefile_Array=(
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/V42Y23-325.gpr
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/V42Y23-325.gpr
    )

manual_alignment_Array=(
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/S13.json
    /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/S15_061343_S13_014988/S15.json
    )

fastq_path_Array=( /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/01_BHK2N3DMXY/outs/fastq_path/H75LKDRX3
                   /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/01_BHK2N3DMXY/outs/fastq_path/H75LKDRX3
                   )

sample_Array=(
    S-13-014988
    S-15-061343
    )

qName=vhighmem
wTime=120:00
cn=10
mem=1000

for i in "${!fastq_path_Array[@]}"; do 
    fastq_path=${fastq_path_Array[$i]}
    image=${image_Array[$i]}
    cyt_image=${cyt_image_Array[$i]}
    slidefile=${slidefile_Array[$i]}
    slide=${slide_Array[$i]}
    area=${area_Array[$i]}
    manual_alignment=${manual_alignment_Array[$i]}
    sample=${sample_Array[$i]}

    JOBFOLDER=${OutDir}/${sample}_with_align
    if [ ! -d $JOBFOLDER ]; then
        mkdir -p $JOBFOLDER
    fi

    JOBNAME=${PIPELINE_NAME}
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
        -n ${cn} \
        -M ${mem} \
        -R rusage[mem=${mem}] \
        -B \
        -N \
        -u xyan4@mdanderson.org \
        /bin/bash -c "
            module load bcl2fastq spaceranger
            spaceranger count --id=${sample} \
                        --transcriptome=/rsrch3/home/genomic_med/xyan4/Software/cellranger/refdata-gex-GRCh38-2020-A \
                        --fastqs=${fastq_path} \
                        --sample=${sample} \
                        --image=${image} \
                        --loupe-alignment=${manual_alignment} \
                        --cytaimage=${cyt_image} \
                        --probe-set=/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/01_data/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                        --slide=${slide} \
                        --slidefile=${slidefile} \
                        --area=${area} "
done

# --loupe-alignment=${manual_alignment} \

