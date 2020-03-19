@echo off
for /f "tokens=1,2 delims=." %%a in ("%1") do (
	set ROOT_FILE=%%a
	set CII_EXTENTION=%%b
)
set EQ=%2
pre_cit %ROOT_FILE%.cii
move %ROOT_FILE%.XSU XSU_MOD
cd XSU_MOD
neutyp_eq %ROOT_FILE%.XSU %EQ% > xsu_data_eq.dat
cd .. 
move XSU_MOD/%ROOT_FILE%.XSU .
cittrs.exe %ROOT_FILE%.cii -SCREEN
copy /Y %ROOT_FILE%.cdb %ROOT_FILE%_eq.cdb
caremdb -opt:export -val:meshflux %ROOT_FILE%_eq.cdb
goto :EOF
