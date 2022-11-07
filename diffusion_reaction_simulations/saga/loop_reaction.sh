#!/bin/bash

# for ITERK in 144
# do
#     for PAT in 176 218 230 199 235 236 105 175 178 183 190 191 228 240 241 227 205 215 002 078 091 127 172 249
#     do
#         echo $PAT "k=" $ITERK
#         # sbatch reaction.slurm $PATID $ITERK gray
#         # sbatch reaction.slurm $PATID $ITERK global
#     done
#     # python -c "t=1; import time; import os; os.system('echo sleep'); os.system('echo ' + str(t)); time.sleep(60*60*t)"
# done
# for PAT in 176 218 230 199 235 236 105 175 178 183 190 191 228 240 241 227 205 215 002 078 091 127 172 249

# those where k = 96, k=144 are already done:
# ['002', '078', '091', '127', '172', '176', '218', '230', '249']
# for PAT in 002 078 091 127 172 176 218 230 249

# those that have to be run again
# for PAT in 091 172 176 190 191 # run without bounds for 10 min time step
# for PAT in 172 176 # run again with time step = 20 min
# for PAT in 091 105 190 191 235 249 # run again with time step = 5 min
# for PAT in 172 175 178 191 190 235 249
# for PAT in 002 078 091 127 172 218 249 # Those where Jplain and Jinit needs to be stored
#199 215 241 # 215 241 # 199 235 236 105 175 178 183 190 191 228 240 241 227 205 215
# for PAT in 178 183 227 235 240

#########################################################
# TODO  1) run for R = 0 for k = 144, 288, 576
# TODO 2) run for R = 5e-5 for k = 96 144, 288, 576
#########################################################

for PAT in 091 176 178
do
    for ITERK in 576
    do
        #########################################################
        # check todo
        #########################################################
        sbatch reaction.slurm $PAT $ITERK
        echo $PAT $ITERK "submitted"
    done
done
