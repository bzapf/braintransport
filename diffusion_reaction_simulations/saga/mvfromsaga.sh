#!/bin/bash

set -e

echo "setting pwd to cd /home/basti/Dropbox\ \(UiO\)/Sleep/"
cd /home/basti/Dropbox\ \(UiO\)/Sleep/


for KF in k144 k288 k96
do
    ITERPATH=iter100
    # for RESNAME in opt_grayalpha_r_tests # opt_r_tests_const # opt_r_tests
    # for RESNAME in opt_r_tests_const
    for RESNAME in opt_r_tests_const_higherBounds
    do
        # for PAT in 002 078 091 127 172 218 249
        for PAT in 172
        # for PAT in 091 191 # 105 172 176 190
        do
            RESPATH=./$PAT/$RESNAME
            mkdir -pv $RESPATH
            SAGADIR=/cluster/projects/nn9279k/bastian/SleepData/$PAT
            FINALPATH=$RESPATH/$ITERPATH/$KF/
            # rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/params $FINALPATH/params

            # echo saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/params
            if [ -d $FINALPATH ]; then
                echo $PAT $FINALPATH "exists, continue"
                continue
            else
                echo "trying to synch" $PAT
                mkdir -pv $RESPATH/$ITERPATH
                mkdir -v $FINALPATH
                { # try to sync
                    echo "try to sync"
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/params $FINALPATH/params
                    echo "synched params"
                } ||
                { # remove path again
                    rm -r $FINALPATH
                    echo $PAT "failed, continue"
                    continue
                    # # synch again to make script exit
                    # rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/params $FINALPATH/params
                }
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/alpha_during_optim.txt $FINALPATH/alpha_during_optim.txt
                    # rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/r_n_during_optim.txt $FINALPATH/r_n_during_optim.txt
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/r_d_during_optim.txt $FINALPATH/r_d_during_optim.txt
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/experimental_data.csv $FINALPATH/experimental_data.csv
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/concs.csv $FINALPATH/concs.csv
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/J_during_optim.txt $FINALPATH/J_during_optim.txt
                    rsync saga:$SAGADIR/$RESNAME/$ITERPATH/$KF/concs_plain.csv $FINALPATH/concs_plain.csv
                    echo "------------------------------------------------------------------------------------"
                    echo $PAT "done"
                    echo "------------------------------------------------------------------------------------"
            fi
            echo $PAT "end of loop"
        done
    done
done