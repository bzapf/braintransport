#!/usr/bin/bash

#SBATCH --job-name=compute_sleep_reduced_ocd
#SBATCH --account=NN9279K
#SBATCH --time=0-04:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --ntasks=1
#SBATCH -o slurmouts/%j.out

## Saga basics:
## 200 standard compute notes with 40 cores and 192 GiB memory each
## 120 standard compute notes with 52 cores and 192 GiB memory each
## 28 medium memory compute nodes, with 40 cores and 384 GiB of memory each
## Nomenclature: core == cpu, task == process

## Run with e.g. --ntasks=1 (16 meshes) --ntasks=1 (32 meshes) or ? (64 meshes)

# Recommended lines from the Sigma2 documentation
# https://documentation.sigma2.no/jobs/job_scripts.html
set -o errexit  
set -o nounset
module --quiet purge # Reset the modules to the system default

module load matplotlib/3.0.0-intel-2018b-Python-3.6.6
module list

# Source common FEniCS installation
# source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf
# Source FEniCS installation with newer scipy
source /cluster/shared/fenics/conf/dolfin-adjoint_for_fenics-2019.1.0.saga.intel.conf

# Create work directory for this job (no need for complicated name,
# will copy results into suitable structure later. The environment
# variable USERWORK is set as /cluster/work/users/meg and the
# SLURM_JOB_ID is given by the job number given from the Slurm system.
workdir=${USERWORK}/modelling-sleep-deprivation-2021/jobs/${SLURM_JOB_ID}
mkdir -pv $workdir

# Set common data-directory
data=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer

# Patient number patient and mesh size n are given as input to job script
patient=$1
n=$2
beta=$3
key=$4

usedti=true
resultfoldername=dtitest
figfolder=FIGURES_CONC_LUT

# usedti=false
# resultfoldername=nodtitest
# figfolder=FIGURES_CONC

# # Copy results back from work to project directory
project_dir=/cluster/projects/nn9279k
results_dir=${project_dir}/meg/modelling-sleep-deprivation-2021/${resultfoldername} # results-red-ocd-b

echo "Job ${SLURM_JOB_ID}: Processing patient ${patient} and mesh size ${n}, beta ${beta} from ${key}."

# Copy meshes and images for this patient into the work directory
cd $data
if [ -e ${data}/${patient}/mesh ] && [ -e ${data}/${patient}/${figfolder} ];
then
    mkdir -pv ${workdir}/${patient}/mesh
    cp -v ${data}/${patient}/mesh/parenchyma*.h5 ${workdir}/${patient}/mesh/

    mkdir -pv ${workdir}/${patient}/${figfolder}
    cp -v ${data}/${patient}/${figfolder}/* ${workdir}/${patient}/${figfolder}/ 
fi

# Copy the script (and other modules used as imports) into the current
# working directory and go there
codedir=/cluster/home/bazapf/nobackup/braintransport/optimal_velocity
SCRIPT_NAME=compute_ocd_reduced.py
cp -v ${codedir}/${SCRIPT_NAME} ${workdir}/
cp -v ${codedir}/compute_ocdr.py ${workdir}/
cp -v ${codedir}/compute_omt.py ${workdir}/
# cp -v ${codedir}/postprocess_omt.py ${workdir}/
cp -v ${codedir}/postprocess_phi.py ${workdir}/
cd ${workdir}

# Submit the script to the job queue
time srun python3 ${SCRIPT_NAME} ${patient} ${n} ${beta} ${key} ${usedti} ${figfolder} > compute_reduced_ocd_pat${patient}_n${n}_beta${beta}_from_${key}.log



if [ ! -d "${results_dir}/${patient}" ];
then
    echo "${results_dir}/${patient} does not exist, creating structures" 
    mkdir -pv ${results_dir}/${patient}/simulations
    mkdir -pv ${results_dir}/${patient}/mesh
    mkdir -pv ${results_dir}/${patient}/${figfolder}
    cp -v ${workdir}/${patient}/mesh/* ${results_dir}/${patient}/mesh/
    cp -v ${workdir}/${patient}/${figfolder}/* ${results_dir}/${patient}/${figfolder}/ 
fi
cp -rv ${workdir}/${patient}/*results* ${results_dir}/${patient}/simulations/
cp -v ${workdir}/*.log ${results_dir}/${patient}/simulations/
