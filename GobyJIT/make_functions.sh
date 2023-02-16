swig -lua -c++ -Iinclude functions.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o functions.so functions_wrap.cxx -lstdc++ -lm -lluajit
