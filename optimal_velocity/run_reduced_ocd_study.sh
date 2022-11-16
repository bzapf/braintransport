#!/usr/bin/bash

patients=(228)

# alphas=(1e-4)
# keys=(0h)
# keys=(6h 24h)
# alphas=(1.0e-03 1.0e-05)

alphas=(1.0e-04 1.0e-03 1.0e-05)
keys=(0h 6h 24h)

## Coarse mesh: n = 16
# for pat in ${patients[@]}; do
#    for alpha in ${alphas[@]}; do
#        for key in ${keys[@]}; do
#            echo "Submitting job:"
#            echo "sbatch --ntasks 1 launch_reduced_ocd_job .sh ${pat} 16 ${alpha} ${key}"
#            sbatch --ntasks 1 launch_reduced_ocd_job.sh ${pat} 16 ${alpha} ${key}
#        done
#    done
# done
# Total time ~10 min for 50 iterations per patient per alpha, peak memory < 1000 MB

# Standard mesh: n = 32
for pat in ${patients[@]}; do
    for alpha in ${alphas[@]}; do
        for key in ${keys[@]}; do
            echo "Submitting job:"
            echo "sbatch --ntasks 1 launch_reduced_ocd_job.sh ${pat} 32 ${alpha} ${key}"
            sbatch --ntasks 1 launch_reduced_ocd_job.sh ${pat} 32 ${alpha} ${key}
        done
    done
done
## Total time ~ N min for 50 iterations per patient per alpha, memory ~
## 2800 MB at 50 iterations, increasing about 30 MB per iteration
