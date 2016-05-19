#!/bin/bash
make -f Makefile.LNX GPU_ARCH=sm_$1 NDEBUG=$2 $3
