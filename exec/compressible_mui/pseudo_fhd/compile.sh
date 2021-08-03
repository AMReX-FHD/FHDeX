#!/bin/bash

echo "** compiling pseudo_fhd.cpp"
mpic++ -std=c++11 -I../../../../MUI pseudo_fhd.cpp -o pseudo_fhd.bin
