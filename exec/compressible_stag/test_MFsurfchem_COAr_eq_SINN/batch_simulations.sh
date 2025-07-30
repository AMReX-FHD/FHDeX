#!/bin/bash

#Run simulations with seeds from 1 to 400
for i in {1..400}
do
    echo "Running seed = $i"
    ../main3d.gnu.MPI.ex inputs_fhd_stag sample_output_flag=0 seed=$i | grep DATA | tee "data${i}.txt"
done

