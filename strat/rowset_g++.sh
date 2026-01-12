#!/bin/bash
g++$3 -Ofast -march=native -DN=$1u $2 rowset.cpp -o rowset_$1.exe
