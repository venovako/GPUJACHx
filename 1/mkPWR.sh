#!/bin/bash
make -f Makefile.PWR GPU_ARCH=$1 NDEBUG=$2 $3
