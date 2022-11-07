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
#SBATCH -o /cluster/home/bazapf/nobackup/sleepCode/scripts/saga/slurmouts/%j.out
# Memory per CPU:
#SBATCH --mem-per-cpu=8G

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

ITERATIONS=100


ARG3=${3:-"gray"}
if [ $ARG3 = "gray" ]
then
    outfolder=opt_grayalpha_r_tests/iter${ITERATIONS}/k${2}
    GLOBALFLAG=--gray
elif [ $ARG3 = "global" ]
then
    outfolder=opt_r_tests_globalalpha/iter${ITERATIONS}/k${2}
    GLOBALFLAG=--global
else
    echo "use either gray or global as 3th argument"
    exit 1
fi



MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
if [ -f $MESHFILE ]; then
echo $MESHFILE "found"
else
echo $MESHFILE "was not found"
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32.h5
echo $MESHFILE "instead"
fi

ARG4=${4:-"nodti"}
if [ $ARG4 = "dti" ]
then
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
echo "arg4 is dti, using mesh with DTI"
DTIFLAG=--dti
else
DTIFLAG=--nodti
fi


ARG5=${5:-"avgt1"}
if [ $ARG5 = "avgt1" ]
then
CONCFOLDER=FIGURES_CONC
else
echo "using FIGUES_CONC_LUT & DTI"
CONCFOLDER=FIGURES_CONC_LUT
outfolder=opt_r_tests_const/iter${ITERATIONS}/k${2}
echo "using FIGURES_CONC_LUT and DTI, setting outfolder to older folder"
echo "outfolder=" $outfolder
MESHFILE=/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/$1/mesh/parenchyma32_with_DTI.h5
DTIFLAG=--dti
fi

OUTFOLDERNAME=opt_r_tests_const_higherBounds


outfolder=$OUTFOLDERNAME/iter${ITERATIONS}/k${2}

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

# Copy script
cp -rv /cluster/home/bazapf/nobackup/sleepCode/scripts/reaction/* ${workdir}/
# cp -rv /cluster/home/bazapf/nobackup/sleepCode/scripts/plots/definitions.py ${workdir}/definitions.py

SCRIPT_NAME=optimal_ra.py
cp -rv /cluster/home/bazapf/nobackup/sleepCode/scripts/simulations/* ${workdir}/

cd ${workdir}
ls ${workdir}

mkdir -p /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/

srun -n 12 python3 ${SCRIPT_NAME} --pat $1 --mesh parenchyma32.h5 \
--sleepdep --exportpath ${workdir} --path_to_files ${workdir} \
--concfolder ${CONCFOLDER} --iter_k $2 --iters $ITERATIONS $GLOBALFLAG $DTIFLAG \
--outfolder ${outfolder} > log_python_srun.txt
 


echo "script done, copying files"
cp -R $USERWORK/$SLURM_JOB_ID/$outfolder/* /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/
cp $USERWORK/$SLURM_JOB_ID/log_python_srun.txt /cluster/projects/nn9279k/bastian/SleepData/$1/$outfolder/log_python_srun.txt

echo 'success'
