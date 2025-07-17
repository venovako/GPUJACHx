#!/bin/bash
#icpx -DNDEBUG -O3 -xHost reverter.cpp -o reverter.exe
#g++ -DNDEBUG -Ofast -march=native reverter.cpp -o reverter.exe
clang++ -DNDEBUG -Ofast -march=native -integrated-as reverter.cpp -o reverter.exe
