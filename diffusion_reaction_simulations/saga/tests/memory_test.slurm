#!/bin/bash

# Job name:
#SBATCH --job-name=memtest
#
# Project:
#SBATCH --account=NN9279K
#
# Wall time limit:
#SBATCH --time 00:20:00
#
### Use standard partition Partition
#
#SBATCH -o ./%j.out
# Memory per CPU:
#SBATCH --mem-per-cpu=1G

# Number of processes:
#SBATCH --ntasks=4

# CPU's per process:
#SBATCH --cpus-per-task=1


## Set up jo environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

module purge
# module load matplotlib/3.0.0-intel-2018b-Python-3.6.6
# module list

# Source common FEniCS installation
# source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf
# Source FEniCS installation with newer scipy
# source /cluster/shared/fenics/conf/dolfin-adjoint_for_fenics-2019.1.0.saga.intel.conf

source /cluster/shared/fenics/conf/fenics-2019.2.0.dev0.saga.intel-2020a-py3.8.conf

# Work directory
# workdir=$USERWORK/$SLURM_JOB_ID
# mkdir -pv $workdir

SCRIPT_NAME=memtestscript.py
# cp -rv ./${SCRIPT_NAME} ${workdir}/

# cd ${workdir}
# ls ${workdir}

srun -n 4 python3 ${SCRIPT_NAME} --test fenics -n $1 > log_$SLURM_JOB_ID.txt
 
echo 'success'
