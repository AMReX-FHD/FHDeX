#!/bin/bash

## Cleanup current directory
clear
unset ALL_PROXY
rm -r plt* stag*

## Select grids
Nsteps=("005" "010" "020")
Spacedim=("064" "128" "256")
dim="2"

# Nsteps=("002" "004" "008")
# Spacedim=("032" "064" "128")
# dim="3"

Visctype=("1" "2" "neg1" "neg2")
Visc=("1" "2" "-1" "-2")

plottype="stagx"
zeros_plt="0000"

## Tool location and name
conv_tool_loc="amrex/Tools/C_util/Convergence"
conv_tool_exec="DiffSameDomainRefinedStag${dim}d.gnu.ex"

for visc_ind in 0 1 2 3
do

    visctype=${Visctype[$visc_ind]}

    # dir_top="Data_Convergence_CrankNicolson/Data_visc_${visctype}/${dim}D"
    dir_top="Data_Convergence_CrankNicolson_GS/Data_visc_${visctype}/${dim}D"

    for grid in 0 1 2
    do
        folder="plots_${Spacedim[$grid]}^${dim}x${Nsteps[$grid]}"
        dir=$dir_top/$folder
        filename=$plottype$zeros_plt${Nsteps[$grid]}
        file[$grid]=$dir/$filename
    done

    ## Run convergence computations
    echo "VISC TYPE = $visctype"
    echo "COARSE COMPUTED ERROR:"
    ~/$conv_tool_loc/$conv_tool_exec reffile=${file[1]} infile1=${file[0]} norm=1
    echo "FINE COMPUTED ERROR:"
    ~/$conv_tool_loc/$conv_tool_exec reffile=${file[2]} infile1=${file[1]} norm=1
    echo " "

done

#END