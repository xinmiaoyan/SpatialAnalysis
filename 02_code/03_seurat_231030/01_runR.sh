#BSUB -J run_R
#BSUB -W 240:00
#BSUB -o /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat
#BSUB -e /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat
#BSUB -cwd /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat
#BSUB -q vlong
#BSUB -u xyan4@mdanderson.org
#BSUB -n 8
#BSUB -M 270
#BSUB -R rusage[mem=270]

cd /rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/02_code/03_seurat_231030

module load R/4.1.0
module load hdf5/1.10.5-hl

Rscript --vanilla 01_S04.R

