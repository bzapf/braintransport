#!/usr/bin/bash


patients=(228)
alphas=(1.0e-03 1.0e-04 1.0e-05)


# Analyze standard mesh results
for pat in ${patients[@]}; do
    for alpha in ${alphas[@]}; do
        echo "Submitting job:"
        echo "sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 6h"
        sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 6h
    done
done

for pat in ${patients[@]}; do
    for alpha in ${alphas[@]}; do
        echo "Submitting job:"
        echo "sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 24h"
        sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 24h
    done
done

for pat in ${patients[@]}; do
    for alpha in ${alphas[@]}; do
        echo "Submitting job:"
        echo "sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 0h"
        sbatch launch_ocd_analysis.sh ${pat} 32 ${alpha} 0h
    done
done

