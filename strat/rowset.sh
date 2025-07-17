#!/bin/bash
#icpx -O3 -xHost -DN=$1u rowset.cpp -o rowset.exe
#g++ -Ofast -march=native -DN=$1u rowset.cpp -o rowset.exe
clang++ -Ofast -march=native -integrated-as -DN=$1u rowset.cpp -o rowset.exe
