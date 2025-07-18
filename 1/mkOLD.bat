@echo off
set GPU_ARCH=%1
shift
set NDEBUG=%1
shift
nmake.exe /NOLOGO /f Makefile.OLD %1 %2
set NDEBUG=
set GPU_ARCH=
