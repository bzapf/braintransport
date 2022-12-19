#!/bin/bash

# Job name:
#SBATCH --job-name=optimal_ra
#
# Project:
#SBATCH --account=NN9279K
#
# Wall time limit:
#SBATCH --time 5-24
#
### Use standard partition Partition
### Use standard partition #SBATCH --partition=bigmem
#
#SBATCH -o /cluster/home/bazapf/nobackup/braintransport/diffusion_reaction_simulations/slurmouts/%j.out
# Memory per CPU:
#SBATCH --mem-per-cpu=4G

# Number of processes:
#SBATCH --ntasks=12

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

ITERATIONS=1


MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
if [ -f $MESHFILE ]; then
echo $MESHFILE "found"
else
echo $MESHFILE "was not found"
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32.h5
echo $MESHFILE "instead"
fi

ARG3=${3:-"nodti"}
if [ $ARG3 = "dti" ]
then
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
echo "arg3 is dti, using mesh with DTI"
DTIFLAG=--dti
else
DTIFLAG=--nodti
fi


ARG4=${4:-"avgt1"}
if [ $ARG4 = "avgt1" ]
then
CONCFOLDER=FIGURES_CONC
else
echo "using FIGUES_CONC_LUT & DTI"
CONCFOLDER=FIGURES_CONC_LUT
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
DTIFLAG=--dti
fi

# OUTFOLDERNAME=PLAINCLEARANCE # 
OUTFOLDERNAME=TOTALCLEARANCE # TEST

outfolder=$OUTFOLDERNAME/k${2}

echo "using CONFOLDER=" $CONCFOLDER

echo "will store to"
echo /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/


# Work directory
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -pv $workdir
mkdir -pv $workdir/$outfolder


mkdir -pv ${workdir}/$1/mesh/

cp -rv $MESHFILE ${workdir}/$1/mesh/parenchyma32.h5


mkdir -pv ${workdir}/$1/${CONCFOLDER}/
cp -rv /cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/${CONCFOLDER}/* ${workdir}/$1/${CONCFOLDER}/

# Copy scripts
cp -rv /cluster/home/bazapf/nobackup/braintransport/diffusion_reaction_simulations/*.py ${workdir}/
# cp -rv /cluster/home/bazapf/nobackup/sleepCode/scripts/plots/definitions.py ${workdir}/definitions.py

SCRIPT_NAME=optimal_ra.py
# cp -rv /cluster/home/bazapf/nobackup/braintransport/diffusion_reaction_simulations/* ${workdir}/

cd ${workdir}
ls ${workdir}

mkdir -p /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/

srun -n 12 python3 -u ${SCRIPT_NAME} --pat $1 --mesh parenchyma32.h5 \
--exportpath ${workdir} --path_to_files ${workdir} \
--concfolder ${CONCFOLDER} --iter_k $2 --iters $ITERATIONS $DTIFLAG --alpha $5 --reaction_rate $6 \
--outfolder ${outfolder} > log_python_srun.txt
 


echo "script done, copying files"
cp -R $USERWORK/$SLURM_JOB_ID/$outfolder/* /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/
cp $USERWORK/$SLURM_JOB_ID/log_python_srun.txt /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/log_python_srun.txt

echo 'success'
