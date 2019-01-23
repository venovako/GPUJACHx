@echo off
for /D %%I in (mmstep,BrentL,colcyc,cycloc,rowcyc,cycwor) do call proc_dir-win.bat %%I
del /F strat.dll
del /F strat.exp
del /F strat.lib
del /F strat.map
link.exe /NOLOGO /VERBOSE /DLL /NOENTRY /MAP /MAPINFO:EXPORTS /RELEASE /OUT:strat.dll *.obj
del /F *.obj
