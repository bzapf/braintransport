#!/bin/bash

# Job name:
#SBATCH --job-name=optimal_ra
#
# Project:
#SBATCH --account=NN9279K
#
# Wall time limit:
#SBATCH --time 0-1
#
### Use standard partition Partition
### Use standard partition #SBATCH --partition=bigmem
#
#SBATCH -o /cluster/home/bazapf/nobackup/sleepCode/scripts/saga/slurmouts/%j.out
# Memory per CPU:
#SBATCH --mem-per-cpu=1G

# Number of processes:
#SBATCH --ntasks=2

# CPU's per process:
#SBATCH --cpus-per-task=1

# ARG5=${5:-"dti"}


## Set up jo environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

module load matplotlib/3.0.0-intel-2018b-Python-3.6.6
module list

# Source common FEniCS installation
# source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf
# Source FEniCS installation with newer scipy
source /cluster/shared/fenics/conf/dolfin-adjoint_for_fenics-2019.1.0.saga.intel.conf

srun -n 2 python3 /cluster/home/bazapf/nobackup/sleepCode/scripts/saga/test_projection.py > testlog.txt
 