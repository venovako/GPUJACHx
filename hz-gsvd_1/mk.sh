#!/bin/bash
make GPU_ARCH=sm_$1 NDEBUG=$2 $3
