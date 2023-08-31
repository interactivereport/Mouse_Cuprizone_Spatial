#!/bin/bash

#Paths
prog_dir="/edgehpc/dept/compbio/projects/SingleCell/programs/dev/spatial_tx/cside"

log_dir=$prog_dir/log/svg_pip

if [[ ! -e $log_dir ]]; then
  mkdir $log_dir
fi


#Job submission

declare -a arr_alt=("cortex" "cortex" "cortex" "cortex" "CC" "CC" "CC" "hippocampus" "hippocampus" "TH")

declare -a arr_ref=("CC" "hippocampus" "TH" "hypoTH" "hippocampus" "TH" "hypoTH" "TH" "hypoTH" "hypoTH")


declare -a grp=("CTL" "CPZ" "RCV")


for j in 0 1 2
do
  for i in {0..9}
  do
  multiline="module load R; \
  Rscript ${prog_dir}/700_cside_pipeline.R \
  --SVG \
  --ave_pop \
  --group ${grp[$j]} \
  --region_ref ${arr_ref[$i]} \
  --region_alt ${arr_alt[$i]} \
  --celltypes en_hpc,en_l5,astrocytes,en_l23,opc,mfol,mol,microglia \
  --log_fc_thresh 0.4 \
  --cell_type_count_thresh 0 \
  --weight_thresh 0 \
  --gene_threshold 0.00005 \
  --max_core 20 \
  --ct_wt /edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result/cell_type_weight_mtx_all.csv \
  --sc_ref /edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/SCRef_1k_ambient_doublet_removed.rds \
  --puck_dir /edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/spacexr_format \
  --qc_spots /edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary/443_ST_QC_scenario5_spots.csv \
  --qc_genes /edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary/443_ST_QC_scenario5_genes.csv \
  --meta_data /edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/metadata.tsv \
  --out_dir /edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/SVG_pip"
  
  srun --job-name=cside -o ${log_dir}/${grp[$j]}_${arr_alt[$i]}_${arr_ref[$i]}_cside_svg_replicates_full_no_filter_QCed.log -p cpu --time='300:00:00' --mem-per-cpu=240G bash -c "$multiline" &
  done
done



