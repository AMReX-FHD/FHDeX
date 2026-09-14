#!/bin/bash

## Cleanup current directory
clear
unset ALL_PROXY
rm -r plt* stag*

# Select grids
Nsteps=("010" "020" "040")
Spacedim=("032" "064" "128")
dims="3"

Visctype=("1" "2" "neg1" "neg2")
Visc=("1" "2" "-1" "-2")

# Choose name of output plotfile to perform convergence tests
plottype="stagx"
zeros_plt="0000"

data_dir="Data_Convergence_CR"
## Tool location and name
conv_tool_loc="../../../amrex/Tools/C_util/Convergence"
conv_tool_exec="DiffSameDomainRefinedStag${dims}d.gnu.ex"

# Loop over viscosity types
for visc_ind in 1
do

    # Create directory name
    visctype=${Visctype[$visc_ind]}
    dir_top="${data_dir}/${dims}D/Data_visc_${visctype}"

    # Loop over grids
    for grid in 0 1 2
    do
        # Collect filename
        folder="plots_${Spacedim[$grid]}^${dims}x${Nsteps[$grid]}"
        dir=$dir_top/$folder
        filename=$plottype$zeros_plt${Nsteps[$grid]}
        file[$grid]=$dir/$filename
    done

    ## Run convergence computations
    echo "VISC TYPE = $visctype"
    echo "COARSE COMPUTED ERROR:"
    $conv_tool_loc/$conv_tool_exec reffile=${file[1]} infile1=${file[0]} norm=1
    echo "FINE COMPUTED ERROR:"
    $conv_tool_loc/$conv_tool_exec reffile=${file[2]} infile1=${file[1]} norm=1
    echo " "

done

#END