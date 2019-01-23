@echo off
icl.exe /nologo /fast /Qcxx-features /DNDEBUG reverter.cpp
del reverter.obj
