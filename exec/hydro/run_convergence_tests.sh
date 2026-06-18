#!/bin/bash

## Cleanup current directory
clear
unset ALL_PROXY
rm -r plt* stag*

# # Select grids
# Nsteps=("010" "020" "040")
# Spacedim=("064" "128" "256")
# Maxgrid=("032" "064" "128")
# Dt=("800" "400" "200")
# dim="2"

# Select grids
Nsteps=("010" "020" "040")
Spacedim=("032" "064" "128")
Maxgrid=("016" "032" "064")
Dt=("400" "200" "100")
dim="3"

nprocs="4"
make -j${nprocs} DIM=${dim}

Visctype=("1" "2" "neg1" "neg2")
Visc=("1" "2" "-1" "-2")

input_file="inputs_dtrmnistic_convergence_${dim}d"
output_dir="Data_Convergence_CR"
mkdir "${output_dir}"

# Reflect test information in compute_convergence_rate.sh
compute_script="compute_convergence_rate.sh"
sed -i "s/Nsteps=.*/Nsteps=(\"${Nsteps[0]}\" \"${Nsteps[1]}\" \"${Nsteps[2]}\")/" ./$compute_script
sed -i "s/Spacedim=.*/Spacedim=(\"${Spacedim[0]}\" \"${Spacedim[1]}\" \"${Spacedim[2]}\")/" ./$compute_script
sed -i "s/dims=.*/dims=\"${dim}\"/" ./$compute_script
sed -i "s/data_dir=.*/data_dir=\"${output_dir}\"/" ./$compute_script

# Fill prob_hi
probhi_placezeros="000" # placezeros for Spacedim scaling
                        # (i.e. if Spacedim=64 & probhi_placezeros=000,
                        # then probhi=64000.)
grid_probhi="2" # Grid index used to compute probhi
if [ $dim = "2" ]
then
    sed -i "s/prob_hi(1:${dim}) = .*/prob_hi(1:${dim}) = ${Spacedim[$grid_probhi]}${probhi_placezeros}. ${Spacedim[$grid_probhi]}${probhi_placezeros}./" ./$input_file
fi
if [ $dim = "3" ]
then
    sed -i "s/prob_hi(1:${dim}) = .*/prob_hi(1:${dim}) = ${Spacedim[$grid_probhi]}${probhi_placezeros}. ${Spacedim[$grid_probhi]}${probhi_placezeros}. ${Spacedim[$grid_probhi]}${probhi_placezeros}./" ./$input_file
fi

# Loop over viscosity types
for visc_ind in 1
do

    visctype=${Visctype[$visc_ind]}

    ## Make directories if not already existing
    mkdir "${output_dir}/${dim}D"
    mkdir "${output_dir}/${dim}D/Data_visc_${visctype}"
    dir_top="${output_dir}/${dim}D/Data_visc_${visctype}"

    sed -i "s/visc_type = .*/visc_type = ${Visc[$visc_ind]}/" ./$input_file

    for grid in 0 1 2
    do

        # Replace variable values in inputs file
        sed -i "s/fixed_dt = .*/fixed_dt = ${Dt[$grid]}/" ./$input_file
        sed -i "s/max_step = .*/max_step = ${Nsteps[$grid]}/" ./$input_file
        sed -i "s/plot_int = .*/plot_int = ${Nsteps[$grid]}/" ./$input_file
        # sed -i "s/plot_int = .*/plot_int = 1/" ./$input_file

        if [ $dim = "2" ]
        then
            sed -i "s/n_cells(1:${dim}) = .*/n_cells(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file
            sed -i "s/max_grid_size(1:${dim}) = .*/max_grid_size(1:${dim}) = ${Maxgrid[$grid]} ${Maxgrid[$grid]}/" ./$input_file
        fi

        if [ $dim = "3" ]
        then
            sed -i "s/n_cells(1:${dim}) = .*/n_cells(1:${dim}) = ${Spacedim[$grid]} ${Spacedim[$grid]} ${Spacedim[$grid]}/" ./$input_file
            sed -i "s/max_grid_size(1:${dim}) = .*/max_grid_size(1:${dim}) = ${Maxgrid[$grid]} ${Maxgrid[$grid]} ${Maxgrid[$grid]}/" ./$input_file
        fi
        ###########################################################################

        folder="plots_${Spacedim[$grid]}^${dim}x${Nsteps[$grid]}"
        dir=$dir_top/$folder

        # Cleanup
        mkdir $dir
        rm -r $dir/stag* $dir/plt*

        # If wish to remove/cleanup previous directory:
        # rm -r $dir

        # Run various inputs files & store plot files in directory
        mpiexec -n ${nprocs} ./main${dim}d.gnu.MPI.ex ${input_file}
        # ./main${dim}d.gnu.MPI.ex ${input_file}

        # Plot result after each run
        # amrvis${dim}d -a plt*

        # Move plot files to specified directory
        mv plt* stag* $dir

    done
done

#END