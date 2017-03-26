#!/bin/bash
./proc_dir-mac.sh mmstep
./proc_dir-mac.sh BrentL
./proc_dir-mac.sh colcyc
./proc_dir-mac.sh cycloc
./proc_dir-mac.sh rowcyc
./proc_dir-mac.sh cycwor
rm -f strat.dylib
rm -f strat.map
ld *.o -dylib -macosx_version_min 10.12 -map strat.map -o strat.dylib
rm -f *.o
