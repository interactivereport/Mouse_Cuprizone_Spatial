#!/bin/bash

#Paths
prog_dir="/edgehpc/dept/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE"
data_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed"
meta_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed"
out_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/RCV_vs_CTRL_ambient_doublet_removed"
log_dir="${out_dir}/log"

mkdir -p $out_dir
mkdir -p $log_dir

#Job submission
for celltype in en_l6b en_l5 en_hpc in_sst en in_pvalb en_thal en_l56 in_vip in_hypo en_l6 in_lamp5 en_l23 en_l25 astrocytes opc nfol endothelial cop microglia mol mfol
do
multiline="module load R; \
Rscript $prog_dir/402_nebula_DE_ambient_doublet_removed.R \
--in_dir $data_dir \
--gene_info_file genes.tsv \
--meta_file $meta_dir/meta_format.tsv \
--count_file matrix.mtx \
--celltype_col Celltype_annotation_clean \
--sample_col Number \
--treatment_col Group \
--cluster $celltype \
--ref_group Control \
--alt_group Recovery \
--lib_size_high 100000 \
--minimum.cells.per.gene.type and \
--minimum.cells.per.subject 10 \
--out_dir $out_dir \
--out_prefix nebula_hl_rcv_vs_ctrl \
--volcano
"
srun --job-name=nebula_hl_${celltype}_rcv_vs_ctrl -o $log_dir/nebula_hl_${celltype}_rcv_vs_ctrl.log -p cpu --time='300:00:00' --mem-per-cpu=512G bash -c "$multiline" &
done