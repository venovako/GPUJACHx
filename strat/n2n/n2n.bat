@echo off
rem Usage: n2n.bat N NNNNN MMMMN BASENAME
echo // %2 > %2.tmp
echo #ifndef C >> %2.tmp
echo #define C 2u >> %2.tmp
echo #endif // !C >> %2.tmp
echo #ifdef N >> %2.tmp
echo #ifndef Sn >> %2.tmp
echo #define Sn ((N) - 1u) >> %2.tmp
echo #endif // !Sn >> %2.tmp
echo #ifndef Pn >> %2.tmp
echo #define Pn ((N) / 2u) >> %2.tmp
echo #endif // !Pn >> %2.tmp
echo #else // !N >> %2.tmp
echo #error N must be defined! >> %2.tmp
echo #endif // N >> %2.tmp
echo #ifndef USE_STRAT_ARRAY_DECLARATOR >> %2.tmp
echo #define USE_STRAT_ARRAY_DECLARATOR >> %2.tmp
echo #endif // !USE_STRAT_ARRAY_DECLARATOR >> %2.tmp
echo #include "%2.h" >> %2.tmp
echo #ifndef genstrat >> %2.tmp
echo #define genstrat(name) %4##name >> %2.tmp
echo #endif // !genstrat >> %2.tmp
echo #ifndef stratN >> %2.tmp
echo #define stratN genstrat(%2) >> %2.tmp
echo #endif // !stratN >> %2.tmp
echo #ifndef stratM >> %2.tmp
echo #define stratM genstrat(%3) >> %2.tmp
echo #endif // !stratM >> %2.tmp
copy /v /y /a %2.tmp /a + %0.cpp %2.cpp /a
icl.exe /nologo /fast /Qcxx-features /DN=%1u %2.cpp
%2.exe %4 > %3.h
del %2.exe
del %2.obj
del %2.cpp
del %2.tmp
