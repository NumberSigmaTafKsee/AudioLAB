swig -lua -c++ vadpwoscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vadpwoscillators.so vadpwoscillators_wrap.cxx -lstdc++ -lm -lluajit
