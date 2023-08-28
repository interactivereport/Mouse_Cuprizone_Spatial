#!/bin/bash

#Directories
prog_dir="/edgehpc/dept/compbio/projects/SingleCell/programs/dev/spatial_tx/spatalk"

#ST data directories
st_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS"
st_qc_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary"

#ST deconvolution directory
c2l_dir="/edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result"

#Single cell reference data directories
sc_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed"

#Single cell downsample reference data directories
sc_ds_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside"

#Manual config file directory
man_input_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/spatalk"

#Output directory
out_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/spatalk/901_spatalk_pipeline_regional_v6_cc_2k_ds_grp_ref"
log_dir=$out_dir/log

#Make log directory; update the permissions of results
mkdir -p $log_dir

#cpz : 106A 072A 074C
#ctrl : 108D 107A 076D
#rcv : 072C 076B 108B

#Sections in Cuprizone
for section in 106A 072A 074C
do
multiline="module load R/4.2.0-foss-2021b; \
Rscript $prog_dir/900_spatalk_pipeline_regional_v6.R \
--st_meta_data $st_dir/metadata.tsv \
--qc_spots $st_qc_dir/443_ST_QC_scenario5_spots.csv \
--qc_genes $st_qc_dir/443_ST_QC_scenario5_genes.csv \
--section $section \
--region_column annotation \
--region_subset CC \
--ct_wt $c2l_dir/cell_type_weight_mtx_all.csv \
--sc_ref_10x_dir $sc_dir \
--sc_celltype_meta_data $sc_dir/meta_format_clean.csv \
--sc_downsample \
--sc_downsample_path $sc_ds_dir/SCRef_2k_ambient_doublet_removed_unassigned_celltypes_excluded_seed_Cuprizone.rds \
--ct_decon \
--n_perm 1000 \
--min_percent 0.1 \
--cci \
--man_sender_receiver $man_input_dir/805_manual_sender_receiver.csv \
--man_lr_path $man_input_dir/805_manual_lrs.csv \
--res_out \
--res_path $out_dir/${section}_c2l_wt_mtx_spatalk_dec_celltype_cci.rds \
--out_dir $out_dir \
--max_core 12 \
--seed 2023
"
srun --job-name=${section}_spatalk -o $log_dir/901_${section}_cc_spatalk_pipeline_v6_regional_cc_2k_ds_seed_group_ref.log -p cpu --time='300:00:00' --mem-per-cpu=264G bash -c "$multiline" &
done
