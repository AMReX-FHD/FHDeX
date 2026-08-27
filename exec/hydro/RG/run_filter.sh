#!/bin/bash

# Default parameter file
PARAM_FILE=${1:-params.txt}

# Read parameters from file
while IFS='=' read -r key value; do
    # Skip empty lines and comments
    [[ -z "$key" || "$key" =~ ^#.*$ ]] && continue
    declare "$key=$value"
done < "$PARAM_FILE"

: "${inputs_file:?Set inputs_file to the original hydro inputs file}"
: "${target_dir:?Set target_dir to the directory containing the checkpoints}"
: "${mpi_ranks:=8}"
: "${plot_filter:=1}"
: "${plot_fourier:=0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC_PATH="${executable:-$SCRIPT_DIR/main2d.gnu.MPI.ex}"

if [[ "$inputs_file" = /* ]]; then
    INPUTS_PATH="$inputs_file"
else
    INPUTS_PATH="$(realpath "$inputs_file")"
fi

if [[ ! -f "$INPUTS_PATH" ]]; then
    echo "Hydro inputs file not found: $INPUTS_PATH" >&2
    exit 1
fi
if [[ ! -x "$EXEC_PATH" ]]; then
    echo "Postprocessor executable not found: $EXEC_PATH" >&2
    exit 1
fi

# Move into the target directory
cd "$target_dir" || { echo "Failed to cd into $target_dir"; exit 1; }

# Print out the available checkpoint files in this directory
echo "========================================="
echo "Available chk directories in $PWD:"
ls -d chk* 2>/dev/null || echo "  No chk directories found!"
echo "========================================="


# Outer loop: iterate mathematically over the checkpoint numbers
for (( chk_num=chk_start; chk_num<=chk_end; chk_num+=chk_inc )); do

    # Hydro checkpoints use a seven-digit step suffix.
    printf -v CHK_NAME 'chk%07d' "$chk_num"
    output_step=$((chk_num + 1))
    echo "========================================="
    echo "Processing checkpoint: $CHK_NAME"
    echo "========================================="

    # Inner loop: logarithmic increment using integer arithmetic
    for (( kmax=kmax_start; kmax<=kmax_end; kmax=(kmax * kmax_inc) / 10 )); do
        FILTER_FILE="filtered_${output_step}_${kmin}_${kmax}"
        FOURIER_FILE="filtered_fourier_${output_step}_${kmin}_${kmax}"
        run_needed=0

        if (( plot_filter != 0 )) && [[ ! -e "$FILTER_FILE" ]]; then
            run_needed=1
        fi
        if (( plot_fourier != 0 )) && [[ ! -e "$FOURIER_FILE" ]]; then
            run_needed=1
        fi
        if (( plot_filter == 0 && plot_fourier == 0 )); then
            run_needed=1
        fi

        if (( run_needed == 0 )); then
            echo "  Requested outputs already exist. Skipping."
        else
            echo "  Running kmin = $kmin and kmax = $kmax on $CHK_NAME"

            mpirun -n "$mpi_ranks" "$EXEC_PATH" "$INPUTS_PATH" \
                restart_file="$CHK_NAME" kmin="$kmin" kmax_list="$kmax" \
                plot_fourier="$plot_fourier" plot_filter="$plot_filter"
        fi
    done
done
