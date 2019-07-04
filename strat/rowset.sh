#!/bin/bash
#icpc -fast -DN=$1u rowset.cpp -o rowset.exe
#gcc -Ofast -march=native -DN=$1u rowset.cpp -o rowset.exe
clang++ -Ofast -march=native -integrated-as -DN=$1u rowset.cpp -o rowset.exe
