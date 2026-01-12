@ECHO OFF
icx.exe /nologo /O3 /QxHost /DN=%1u %2 /Ferowset_%1.exe rowset.cpp /link /RELEASE
