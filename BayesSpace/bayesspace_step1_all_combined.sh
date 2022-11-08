#!/bin/bash
#SBATCH --job-name=bayesspace
#SBATCH --array=1
#SBATCH --time=60:00:00
#SBATCH --mem=100000
#SBATCH --mail-user=spiya@biogen.com
#SBATCH --mail-type=END,FAIL

module use /usr/prog/modules/all
module load R/4.1.3-foss-2021b

cd /edgehpc/dept/compbio/projects/TST11523/mouse_cpz

R CMD BATCH --no-save --no-restore ./bayesspace_step1_all_combined.R bayesspace_step1_${SLURM_ARRAY_TASK_ID}.Rout
