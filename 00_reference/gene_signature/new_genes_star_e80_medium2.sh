#BSUB -J p1_sr
#BSUB -W 23:00
#BSUB -o /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan_HER2/a1/output.txt
#BSUB -e /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan_HER2/a1/error.txt
#BSUB -cwd /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan_HER2/a1/
#BSUB -q e80medium
#BSUB -n 40
#BSUB -M 900
#BSUB -R rusage[mem=900]
#BSUB -B
#BSUB -N
#BSUB -u kcho2@mdanderson.org

module load python
module load hdf5/1.10.5-hl

source activate
conda activate STAR_beta
                  
prefix="/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan_HER2/a1/"
device="cpu"  # "cuda" or "cpu"
#n_genes=3000

python /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/Region_prediction_by_imputation.py ${prefix} /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/TRAV.txt ${prefix}marker_plot_TCR/TRAV/ ${prefix}signature_TCR/TRAV.png
