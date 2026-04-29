#!/bin/bash
#
#SBATCH -J hsp90-d13b-d15b
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --no-requeue
#SBATCH -t 55:00:00

jobname=hsp90-d13b-d15b

. ~/miniconda3/bin/activate atm8.1
echo "Running on $(hostname)"

if [ ! -f ${jobname}_0.xml ]; then
   python /home/users/jwu/bin/atm8.1/AToM-OpenMM/rbfe_structprep.py ${jobname}_asyncre.cntl || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
python /home/users/jwu/bin/atm8.1/AToM-OpenMM/rbfe_explicit.py ${jobname}_asyncre.cntl
