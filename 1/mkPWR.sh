#!/bin/bash
make -f Makefile.PWR GPU_ARCH=sm_$1 NDEBUG=$2 $3
