# CSIDE analysis codes

This folder contains codes used to run the spacexr package "CSIDE" DE analysis on Biogen mouse cuprizone data.

## Running codes

* External inputs (generated client-side)

  * /edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result/cell_type_weight_mtx_all.csv

* Codes to generate internal inputs (assumes ambient RNA and doublet removed data is available)

  * 100-1_make_snrnaseq_reference_1k_ambient_doublet_removed_unassigned_celltypes_excluded.R
  
  _note: not seeded, input file required_: 
  /edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/SCRef_1k_ambient_doublet_removed_unassigned_celltypes_excluded.rds/
  SCRef_1k_ambient_doublet_removed_unassigned_celltypes_excluded.rds
  
  * 300_CSIDE_make_data.R
  
  * 443_ST_QC_summary.R
 
* Codes to run CSIDE analysis

  * 700_cside_pipeline.R
 
* Example HPC submission script

  * 701_srun_cside_pipeline.sh