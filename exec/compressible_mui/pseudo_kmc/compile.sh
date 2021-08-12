#!/bin/bash

echo "** compiling pseudo_kmc.cpp"
mpic++ -std=c++11 -I../../../../MUI pseudo_kmc.cpp -o pseudo_kmc.bin
