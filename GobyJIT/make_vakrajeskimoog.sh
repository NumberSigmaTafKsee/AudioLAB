swig -lua -c++ vakrajeskimoog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vakrajeskimoog.so vakrajeskimoog_wrap.cxx -lstdc++ -lm -lluajit
