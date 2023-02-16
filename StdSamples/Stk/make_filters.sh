swig -lua -c++ -Iinclude Filters.i
gcc -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Filters.so Filters_wrap.cxx  \
-lstdc++ -lm -lluajit
