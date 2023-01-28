swig -lua -c++ -Iinclude audiodsp_functions.i
gcc -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Functions.so audiodsp_functions_wrap.cxx  \
-lstdc++ -lm -lluajit
