#!/bin/bash
make -f Makefile.MAC GPU_ARCH=$1 NDEBUG=$2 $3
