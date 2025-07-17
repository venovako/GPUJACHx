@echo off
rem Usage: asymptoter.bat N ????? BASENAME IsCyclic
rem E.g., ????? = 00032 when N = 32
rem IsCyclic = 1u for a cyclic, or 0u for a quasi-cyclic strategy
echo // %2 > %3%2.tmp
echo #ifdef N >> %3%2.tmp
echo #ifndef S >> %3%2.tmp
echo #define S ((N) - %4) >> %3%2.tmp
echo #endif // !S >> %3%2.tmp
echo #ifndef P >> %3%2.tmp
echo #define P ((N) / 2u) >> %3%2.tmp
echo #endif // !P >> %3%2.tmp
echo #ifndef C >> %3%2.tmp
echo #define C 2u >> %3%2.tmp
echo #endif // !C >> %3%2.tmp
echo #else // !N >> %3%2.tmp
echo #error N must be defined! >> %3%2.tmp
echo #endif // N >> %3%2.tmp
echo #ifndef USE_STRAT_ARRAY_DECLARATOR >> %3%2.tmp
echo #define USE_STRAT_ARRAY_DECLARATOR >> %3%2.tmp
echo #endif // !USE_STRAT_ARRAY_DECLARATOR >> %3%2.tmp
echo #include "../%3/%2.h" >> %3%2.tmp
echo #ifndef strat >> %3%2.tmp
echo #define strat %3%2 >> %3%2.tmp
echo #else // strat >> %3%2.tmp
echo #error strat already defined >> %3%2.tmp
echo #endif // !strat >> %3%2.tmp
copy /v /y /a %3%2.tmp /a + %0.cpp %3%2.cpp /a
icx.exe /nologo /O3 /QxHost /DN=%1u %3%2.cpp /link /RELEASE
%3%2.exe %3%2
del %3%2.exe
del %3%2.obj
del %3%2.cpp
del %3%2.tmp
del texput.aux
del texput.log
pdftk.exe %3%2_??.pdf cat output %3%2.pdf verbose dont_ask
