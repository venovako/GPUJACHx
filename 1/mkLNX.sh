#!/bin/bash
GCCMAJORVER=`gcc -dumpversion | cut -d. -f1`
if [ "$GCCMAJORVER" != "4" ]
then
	if [ "$GCCMAJORVER" != "9" ]
	then
		make -f Makefile.LNX GPU_ARCH=$1 NDEBUG=$2 STD=14 $3
	else
		make -f Makefile.LNX GPU_ARCH=$1 NDEBUG=$2 STD=17 $3
	fi
	
else
	make -f Makefile.LNX GPU_ARCH=$1 NDEBUG=$2 $3
fi
unset GCCMAJORVER
