% http://tdm-gcc.tdragon.net/download 
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');

% compile source code for Simulink
mex -g -v -IC:\Users\Yashren\work\orbit-experiments\sgp4 sgp4\aSGP4.cpp sgp4\asgp4_util.cpp sgp4\sgp4io.cpp sgp4\sgp4ext.cpp sgp4\sgp4unit.cpp