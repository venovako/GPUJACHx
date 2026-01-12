#!/bin/bash
clang++ -Ofast -march=native -integrated-as -DN=$1u $2 rowset.cpp -o rowset_$1.exe
