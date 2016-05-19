#!/bin/bash
./proc_dir-lnx.sh mmstep
./proc_dir-lnx.sh BrentL
./proc_dir-lnx.sh colcyc
./proc_dir-lnx.sh cycloc
./proc_dir-lnx.sh rowcyc
./proc_dir-lnx.sh cycwor
rm -f strat.so
rm -f strat.map
ld -shared -Map=strat.map -o strat.so *.o
rm -f *.o
