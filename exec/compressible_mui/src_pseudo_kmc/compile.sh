#!/bin/bash

echo "** compiling receiver.cpp"
mpic++ -std=c++11 -I../../../../MUI pseudo_kmc.cpp -o ../pseudo_kmc.bin
