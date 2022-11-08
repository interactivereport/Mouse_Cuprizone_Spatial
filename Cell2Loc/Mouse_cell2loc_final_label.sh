#!/bin/bash

#SBATCH --job-name=Mouse_cell2loc
#SBATCH --array=1-28
#SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=300000
##SBATCH --exclusive
#SBATCH --mail-user=sarbottam.piya@biogen.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o ./error_report/slurm-%A_%a.out

module load anaconda3

conda activate /home/spiya/.conda/envs/cell2loc_env
cd /edgehpc/dept/compbio/projects/TST11523/mouse_cpz


sid=$SLURM_ARRAY_TASK_ID
mouse_id=$(awk -v rowno=$sid 'NR == rowno {print $1}' ID.map.txt)

jupyter nbconvert --to notebook --execute Mouse_cell2loc_final_label.ipynb --output=Mouse_cell2loc_final_label_${mouse_id}.ipynb
