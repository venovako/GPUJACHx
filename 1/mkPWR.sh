#!/bin/bash
GCCMAJORVER=`gcc -dumpversion | cut -d. -f1`
if [ "$GCCMAJORVER" != "4" ]
then
	make -f Makefile.PWR GPU_ARCH=$1 NDEBUG=$2 STD=14 $3
else
	make -f Makefile.PWR GPU_ARCH=$1 NDEBUG=$2 $3
fi
unset GCCMAJORVER
