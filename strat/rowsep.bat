@ECHO OFF
icx.exe /nologo /Qopenmp /O3 /QxHost /DN=%1u %2 /Ferowsep_%1.exe rowsep.cpp /link /RELEASE
