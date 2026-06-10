#!/bin/bash

Nx=512
Ny=512
SEED=1
MAX_GRID_SIZE=64

FIXED_DT=1.9e-4
MAX_STEP=100000
PLOT_INT=500

nu=1.95e-3
CELL_DEPTH=100000000.


RAW_DIR="snapshots_Nx${Nx}_Ny${Ny}_seed${SEED}_mgs${MAX_GRID_SIZE}_dt${FIXED_DT}_nu${NU}_depth${CELL_DEPTH}"

RUN_DIR="${RAW_DIR//./p}"

PLT_DIR="${RUN_DIR}/plt"


mpirun -n 16 ./main2d.gnu.MPI.ex inputs_kh_2d plot_base_name="${PLT_DIR}" n_cells="${Nx} ${Ny}" seed="${SEED}" max_grid_size="${MAX_GRID_SIZE} ${MAX_GRID_SIZE}" fixed_dt="${FIXED_DT}" max_step="${MAX_STEP}" plot_int="${PLOT_INT}" visc_coef="${nu}" cell_depth="${CELL_DEPTH}"