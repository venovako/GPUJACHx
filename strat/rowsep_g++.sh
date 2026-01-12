#!/bin/bash
g++$3 -fopenmp -Ofast -march=native -DN=$1u $2 rowsep.cpp -o rowsep_$1.exe
