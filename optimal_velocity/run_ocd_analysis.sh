#!/usr/bin/bash

# patients=(199 235 236 105 175 178 183 190 191 228 240 241 227 205 215)

# patients=(105)

# patients=(199 235 236 175 178 183) # 
patients=(190 191 228 240 241 227 205 215)
alphas=(1.0e-03 1.0e-04 1.0e-05)

patients=(190 191 228 240 241 227 205 215)
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

