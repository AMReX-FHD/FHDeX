#!/bin/bash

NRUN=64

for ((i=1;i<=$NRUN;i++))
do
    echo RUN$i
    cd RUN$i
    cp ../havg.py .
    python havg.py
    cp ../mass_cons.py .
    python mass_cons.py

    sed -i ' 1 s/.*/#&/' res.mass_spec1_final_vert
    sed -i ' 1 s/.*/#&/' res.mass_spec2_final_vert
    sed -i ' 1 s/.*/#&/' res.mass_spec3_final_vert
    sed -i ' 1 s/.*/#&/' res.mass_spec4_final_vert

    cd ..
done
