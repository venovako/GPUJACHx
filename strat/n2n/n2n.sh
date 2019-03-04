#!/bin/bash
#Usage: ./n2n.sh N NNNNN MMMMM BASENAME
echo "// $2" > $2.tmp
echo "#ifndef C" >> $2.tmp
echo "#define C 2u" >> $2.tmp
echo "#endif // !C" >> $2.tmp
echo "#ifdef N" >> $2.tmp
echo "#ifndef Sn" >> $2.tmp
echo "#define Sn ((N) - 1u)" >> $2.tmp
echo "#endif // !Sn" >> $2.tmp
echo "#ifndef Pn" >> $2.tmp
echo "#define Pn ((N) / 2u)" >> $2.tmp
echo "#endif // !Pn" >> $2.tmp
echo "#else // !N" >> $2.tmp
echo "#error N must be defined!" >> $2.tmp
echo "#endif // N" >> $2.tmp
echo "#ifndef USE_STRAT_ARRAY_DECLARATOR" >> $2.tmp
echo "#define USE_STRAT_ARRAY_DECLARATOR" >> $2.tmp
echo "#endif // !USE_STRAT_ARRAY_DECLARATOR" >> $2.tmp
echo "#include \"$2.h\"" >> $2.tmp
echo "#ifndef genstrat" >> $2.tmp
echo "#define genstrat(name) $4##name" >> $2.tmp
echo "#endif // !genstrat" >> $2.tmp
echo "#ifndef stratN" >> $2.tmp
echo "#define stratN genstrat($2)" >> $2.tmp
echo "#endif // !stratN" >> $2.tmp
echo "#ifndef stratM" >> $2.tmp
echo "#define stratM genstrat($3)" >> $2.tmp
echo "#endif // !stratM" >> $2.tmp
cat $2.tmp $0.cpp > $2.cpp
g++ -Ofast -march=native -DN=$1u $2.cpp -o $2.exe -lm
./$2.exe $4 > $3.h
unix2dos $3.h
rm $2.exe
rm $2.cpp
rm $2.tmp