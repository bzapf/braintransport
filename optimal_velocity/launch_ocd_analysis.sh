#!/usr/bin/bash

#SBATCH --job-name=sleep_red_ocd_analysis
#SBATCH --account=NN9279K
#SBATCH --time=0-00:30:00
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks=1
#SBATCH -o slurmouts/%j.out
#SBATCH -e slurmouts/%j.err

## Saga basics:
## 200 standard compute notes with 40 cores and 192 GiB memory each
## 120 standard compute notes with 52 cores and 192 GiB memory each
## 28 medium memory compute nodes, with 40 cores and 384 GiB of memory each
## Nomenclature: core == cpu, task == process

## Run with  --ntasks=1

# Recommended lines from the Sigma2 documentation
# https://documentation.sigma2.no/jobs/job_scripts.html
set -o errexit  
set -o nounset
module --quiet purge # Reset the modules to the system default

module load matplotlib/3.0.0-intel-2018b-Python-3.6.6
module list

# Source common FEniCS installation
source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

# Create work directory for this job (no need for complicated name,
# will copy results into suitable structure later. The environment
# variable USERWORK is set as /cluster/work/users/meg and the
# SLURM_JOB_ID is given by the job number given from the Slurm system.
workdir=${USERWORK}/modelling-sleep-deprivation-2021/jobs/${SLURM_JOB_ID}
mkdir -pv $workdir
# # Copy results and log back from work to project directory
project_dir=/cluster/projects/nn9279k

# Set common data-directory
resultfoldername=dtitest
data=${project_dir}/meg/modelling-sleep-deprivation-2021/${resultfoldername}

# Patient number patient, mesh size n, beta and (time) key are given as input to script
patient=$1
n=$2
beta=$3
key=$4
echo "Job ${SLURM_JOB_ID}: Processing patient ${patient} and mesh size ${n}, beta ${beta} from ${key}"

# Copy relevant data into the work directory with same structure
mkdir -pv ${workdir}/hdf5
cp -v ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/hdf5/phi.h5 ${workdir}/hdf5/
cp -v ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/hdf5/parenchyma*.h5 ${workdir}/hdf5/
mkdir -pv ${workdir}/pvd

# Copy the script into the current working directory and go there
codedir=/cluster/home/bazapf/nobackup/braintransport/optimal_velocity
SCRIPT_NAME=postprocess_phi.py
cp -v ${codedir}/${SCRIPT_NAME} ${workdir}/
cd ${workdir}

# Submit the script to the job queue
time python3 $SCRIPT_NAME . ${n} > postprocess_phi_pat${patient}_n${n}_beta${beta}_at${key}.log


# cp -rv ${workdir}/hdf5/* ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/hdf5/
# cp -v ${workdir}/pvd/* ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/pvd/
cp -v ${workdir}/*.log ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/
cp -v ${workdir}/*.csv ${data}/${patient}/simulations/results_red_ocd_pat${patient}_n${n}_beta${beta}_${key}/

# Sample run:
# sbatch launch_ocd_analysis.sh 002 16 1.0e-04 6h
