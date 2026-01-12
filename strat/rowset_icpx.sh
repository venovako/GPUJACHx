#!/bin/bash
icpx -O3 -xHost -DN=$1u $2 rowset.cpp -o rowset_$1.exe
