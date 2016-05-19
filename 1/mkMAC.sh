#!/bin/bash
make -f Makefile.MAC GPU_ARCH=sm_$1 NDEBUG=$2 $3
