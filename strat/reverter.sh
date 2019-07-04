#!/bin/bash
#icpc -DNDEBUG -fast reverter.cpp -o reverter.exe
#gcc -DNDEBUG -Ofast reverter.cpp -o reverter.exe
clang++ -DNDEBUG -Ofast -march=native -integrated-as reverter.cpp -o reverter.exe
