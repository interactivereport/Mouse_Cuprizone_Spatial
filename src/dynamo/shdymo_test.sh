#!/bin/bash 
#SBATCH --job-name=dyn
#SBATCH -p gpu
##################SBATCH --gres=gpu:1
#SBATCH -t 50-00:00:00 
#SBATCH --nodes=1
#SBATCH --mem=1T
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=32
################SBATCH -a 0-11
#SBATCH -o dymo.out
#SBATCH -e dymo.err

module load anaconda3 
conda activate /edgehpc/dept/compbio/users/whu1/envs/dynamo_
python dynamo_test_tmp.py
conda deactivate 