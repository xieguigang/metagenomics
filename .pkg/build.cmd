@echo off

SET drive=%~d0
SET R_HOME=%drive%/GCModeller\src\R-sharp\App\net5.0
SET pkg=./metagenomics.zip

%R_HOME%/Rscript.exe --build /src ../ /save %pkg%
%R_HOME%/R#.exe --install.packages %pkg%

pause