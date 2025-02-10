#BSUB -J p1_sr
#BSUB -W 23:00
#BSUB -o /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/output_copy.txt
#BSUB -e /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/error_copy.txt
#BSUB -cwd /rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/
#BSUB -q e80medium
#BSUB -n 40
#BSUB -M 900
#BSUB -R rusage[mem=900]
#BSUB -B
#BSUB -N
#BSUB -u kcho2@mdanderson.org

while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/b_t_den_remain/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/b_t_den_remain.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/caf/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/caf.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/ccl_cxcl/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/ccl_cxcl.txt"


while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/cd4_c1_tregs/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/cd4_c1_tregs.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/cd4_cd8_c4_tstr/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/cd4_cd8_c4_tstr.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/cd8_c1_tex/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/cd8_c1_tex.txt"

while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/cxcl/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/cxcl.txt"


while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/MAPK/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/MAPK.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/proliferating_tumor_cell/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/proliferating_tumor_cell.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/tls_imprint/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/tls_imprint.txt"


while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/tls_ours/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/tls_ours.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/TRAV/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/TRAV.txt"



while read gene; do
gene=$(echo "${gene}" | tr -d '[:space:]')   # remove leading/trailing whitespaces
if [ -e "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" ]; then 
cp "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/cnts-super-plots/${gene}.png" "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/st-proj-beta-master_jan/B1_3/marker_plot/tumor_hypoxia/cnts-super-plots/"
echo "File '${gene}.png' has been copied."
else
  echo "File '${gene}.png' does not exist in the source directory."
fi
done < "/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_gene_list/tumor_hypoxia.txt"

