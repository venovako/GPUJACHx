#!/bin/bash
icpx -qopenmp -O3 -xHost -DN=$1u $2 rowsep.cpp -o rowsep_$1.exe
