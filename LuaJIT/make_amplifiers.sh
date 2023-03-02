swig -lua -c++ amplifiers.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o amplifiers.so amplifiers_wrap.cxx -lstdc++ -lm -lluajit
