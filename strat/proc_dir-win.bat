@echo off
cd %1
echo #include "../hdr_win.h" > %1.c
for %%I in (*.h) do ( echo EXPORT_VAR >> %1.c && echo #include "%%I" >> %1.c )
del /F %1.obj
cl.exe /nologo /Ox /favor:INTEL64 /fp:precise /GF /GS- /GT /MD /bigobj /c %1.c
del /F %1.c
move /Y %1.obj ..
cd ..
