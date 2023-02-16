swig -lua -c++ voltage.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o voltage.so voltage_wrap.cxx -lstdc++ -lm -lluajit
