#!/bin/bash
# E.g., ANIMATE=1 HDF5_FULL="mpi,z"
make GPU_ARCH=sm_$1 NDEBUG=$2 HOST_FLAGS="-ffp-contract=on,-integrated-as" CVG=1 $3
