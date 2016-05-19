@echo off
set GPU_ARCH=sm_%1
shift
set NDEBUG=%1
shift
nmake.exe /NOLOGO /f Makefile.WIN %1 %2
set NDEBUG=
set GPU_ARCH=
