@echo off
icl.exe /nologo /fast /Qcxx-features /DN=%1u rowset.cpp /link /RELEASE
