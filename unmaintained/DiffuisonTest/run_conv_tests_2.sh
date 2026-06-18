#!/bin/bash

## Cleanup current directory
clear
unset ALL_PROXY
rm -r plt* stag*

## Select grids
Nsteps=("005" "010" "020")
Spacedim=("064" "128" "256")
# Dt=("1.0e-3" "0.5e-3" "0.25e-3")
Dt=("20.0e-5" "10.0e-5" "5.0e-5")
dim="2"

# Nsteps=("002" "004" "008")
# Spacedim=("032" "064" "128")
# Dt=("4.0e-4" "2.0e-4" "1.0e-4")
# dim="3"

Visctype=("1" "2" "neg1" "neg2")
Visc=("1" "2" "-1" "-2")

input_file="inputs_${dim}d"

for visc_ind in 0 1 2 3
do

    visctype=${Visctype[$visc_ind]}

    # Make directories if not already existing
    # mkdir "Data_Convergence_CrankNicolson"
    # mkdir "Data_Convergence_CrankNicolson/Data_visc_${visctype}"
    # mkdir "Data_Convergence_CrankNicolson/Data_visc_${visctype}/${dim}D"
    # dir_top="Data_Convergence_CrankNicolson/Data_visc_${visctype}/${dim}D"

    mkdir "Data_Convergence_CrankNicolson_GS"
    mkdir "Data_Convergence_CrankNicolson_GS/Data_visc_${visctype}"
    mkdir "Data_Convergence_CrankNicolson_GS/Data_visc_${visctype}/${dim}D"
    dir_top="Data_Convergence_CrankNicolson_GS/Data_visc_${visctype}/${dim}D"

    sed -i "s/visc_type = .*/visc_type = ${Visc[$visc_ind]}/" ./$input_file

    for grid in 0 1 2
    do

        # Replace variable values in inputs file
        sed -i "s/fixed_dt = .*/fixed_dt = ${Dt[$grid]}/" ./$input_file
        sed -i "s/max_step = .*/max_step = ${Nsteps[$grid]}/" ./$input_file
        sed -i "s/plot_int = .*/plot_int = ${Nsteps[$grid]}/" ./$input_file
        # sed -i "s/plot_int = .*/plot_int = 1/" ./$input_file

        # FIXME: need if statement for higher dimensions
        ###########################################################################
        sed -i "s/n_cells(1:${dim}) = .*/n_cells(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file
        sed -i "s/max_grid_size(1:${dim}) = .*/max_grid_size(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file

        # sed -i "s/n_cells(1:${dim}) = .*/n_cells(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file
        # sed -i "s/max_grid_size(1:${dim}) = .*/max_grid_size(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file
        ###########################################################################

        folder="plots_${Spacedim[$grid]}^${dim}x${Nsteps[$grid]}"
        dir=$dir_top/$folder

        ## Cleanup
        mkdir $dir
        rm -r $dir/stag* $dir/plt*

        ## If wish to remove/cleanup previous directory:
        # rm -r $dir

        ## Run various inputs files & store plot files in directory
        ./main${dim}d.gnu.MPI.ex inputs_${dim}d
        # amrvis${dim}d -a plt*
        mv plt* stag* $dir

    done
done

#END