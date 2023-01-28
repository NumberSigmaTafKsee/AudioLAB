swig -lua -c++ vabandlimitedoscillators.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vabandlimitedoscillators.so vabandlimitedoscillators_wrap.cxx -lstdc++ -lm -lluajit
