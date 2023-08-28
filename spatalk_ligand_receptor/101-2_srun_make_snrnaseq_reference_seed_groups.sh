#!/bin/bash

#Paths
prog_dir="/edgehpc/dept/compbio/projects/SingleCell/programs/dev/spatial_tx/cside/"

log_dir=$prog_dir/log

#Job submission
for grp in "Control" "Cuprizone"  "Recovery"
do
  multiline="module load R; Rscript $prog_dir/100-2_make_snrnaseq_reference_2k_ambient_doublet_removed_unassigned_celltypes_excluded_seed_groups.R $grp"
  srun --job-name=mkref_2k_seed_grps -o $log_dir/mkref_2k_nounassigned_seed_group_$grp.log -p cpu --time='300:00:00' --mem-per-cpu=160G bash -c "$multiline" &
done

