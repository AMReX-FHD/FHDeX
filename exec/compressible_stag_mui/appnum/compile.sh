#!/bin/bash

mpicxx main.cpp -o prog1.ex
mpicxx -DSLURM main.cpp -o prog2.ex
