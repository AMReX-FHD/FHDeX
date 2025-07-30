#!/bin/bash

CHUNK=50000
MAX_STEP=100000
EXEC="../../main3d.gnu.MPI.ex"
INPUT="../inputs_fhd_stag"

for i in {1..400}; do
    echo ">>> Starting trajectory $i"

    TRAJ_DIR="traj_$i"
    mkdir -p "$TRAJ_DIR"
    cd "$TRAJ_DIR" || exit 1

    STEP=0
    > "../data${i}.txt"

    while [ $STEP -lt $MAX_STEP ]; do
        echo "  Running step $STEP to $((STEP + CHUNK))"

        if [ $STEP -eq 0 ]; then
            $EXEC $INPUT seed=$i restart=-1 max_step=$CHUNK chk_int=$CHUNK sample_output_flag=0 | grep DATA | tee -a "../data${i}.txt"
        else
            CHK_NAME=$(printf "chk%09d" $STEP)
            if [ ! -d "$CHK_NAME" ]; then
                echo "Checkpoint $CHK_NAME missing! Aborting trajectory $i."
                break
            fi
            $EXEC $INPUT restart=$STEP amr.restart_file=$CHK_NAME max_step=$((STEP + CHUNK)) chk_int=$CHUNK sample_output_flag=0 | grep DATA | tee -a "../data${i}.txt"
        fi

        STEP=$((STEP + CHUNK))
    done

    cd ..
done



