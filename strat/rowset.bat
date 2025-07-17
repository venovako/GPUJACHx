@echo off
icpx.exe /nologo /O3 /QxHost /DN=%1u rowset.cpp /link /RELEASE
