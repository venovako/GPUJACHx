#!/bin/bash
make -f Makefile.LNX GPU_ARCH=$1 NDEBUG=$2 $3
