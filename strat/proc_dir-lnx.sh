#!/bin/bash
pushd $1
echo \#include \"../hdr_unx.h\" > $1.c
for I in *.h
do
    echo EXPORT_VAR >> $1.c
    echo \#include \"$I\" >> $1.c
done
rm -f $1.o
gcc -O3 -march=native -fPIC -c $1.c
rm -f $1.c
mv $1.o ..
popd
