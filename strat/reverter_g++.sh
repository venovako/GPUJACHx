#!/bin/bash
g++$1 -DNDEBUG -Ofast -march=native reverter.cpp -o reverter.exe
