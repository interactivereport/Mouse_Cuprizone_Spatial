# SpaTalk analysis codes

This folder contains codes used to run SpaTalk analysis on Biogen mouse cuprizone data.

## Running codes

* External inputs (generated client-side)

  * /edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result/cell_type_weight_mtx_all.csv

* Codes to generate internal inputs (assumes ambient RNA and doublet removed data is available)

  * 710_mouse_cuprizone_ambient_doublet_removed_meta_format.R

  * 443_ST_QC_summary.R (located in Mouse_Cuprizone_Spatial/cside/ set of codes)

  * 100-2_make_snrnaseq_reference_2k_ambient_doublet_removed_unassigned_celltypes_excluded_seed_groups.R

  * 805_manual_spatalk_config_files.R
 
* Codes to run SpaTalk analysis

  * 000_spatalk_aux.R
 
  * 900_spatalk_pipeline_regional_v6.R
 
* Example HPC submission script

  * 901_srun_spatalk_pipeline_regional_v6_cc_2k_ds_seed_group_ref_cpz.sh