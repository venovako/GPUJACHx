#!/bin/bash
make GPU_ARCH=sm_$1 NDEBUG=$2 HDF5=/usr/local HDF5_FULL="mpi,z" HOST_FLAGS="-ffp-contract=on,-integrated-as" CVG=1 $3
