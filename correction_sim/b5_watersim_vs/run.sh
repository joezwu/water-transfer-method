#!/bin/bash
#
#SBATCH -J b5
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --no-requeue
#SBATCH -t 55:00:00

jobname=b5

. /home/users/jwu/miniconda3/bin/activate atm8.1
echo "Running on $(hostname)"

python /home/users/jwu/bin/atm8.1/AToM-OpenMM/md_ligand_water.py ${jobname}
