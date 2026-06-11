#!/bin/bash

source config.sh

RAW_DIR="snapshots_Nx${Nx}_Ny${Ny}_seed${SEED}_mgs${MAX_GRID_SIZE}_dt${FIXED_DT}_visc_coef${VISC_COEF}_noise_off_step${NOISE_OFF_STEP}_depth${CELL_DEPTH}_k_B${K_B}_T_init${T_INIT}"

RUN_DIR="${RAW_DIR//./p}"

PLT_DIR="${RUN_DIR}/plt"


mpirun -n 32 ./main2d.gnu.MPI.ex inputs_kh_2d plot_base_name="${PLT_DIR}" n_cells="${Nx} ${Ny}" seed="${SEED}" max_grid_size="${MAX_GRID_SIZE} ${MAX_GRID_SIZE}" fixed_dt="${FIXED_DT}" max_step="${MAX_STEP}" plot_int="${PLOT_INT}" visc_coef="${VISC_COEF}" cell_depth="${CELL_DEPTH}" noise_off_step="${NOISE_OFF_STEP}" k_B="${K_B}" T_init="${T_INIT}"