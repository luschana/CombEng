
%mex -setup:C:\Users\luschana\AppData\Roaming\MathWorks\MATLAB\R2016a\mex_C_win64.xml C
%mex -setup:'C:\Program Files\MATLAB\R2016a\bin\win64\mexopts\mingw64_g++.xml' C++

mex -v -g mexengine.cpp combustionengine.cpp cylinder.cpp definitions.cpp ecu.cpp environment.cpp gascomponent.cpp injector.cpp oil.cpp shomate.cpp valve.cpp
