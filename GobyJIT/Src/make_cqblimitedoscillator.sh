swig -lua -c++ -Iinclude cqblimitedoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cqblimitedoscillator.so cqblimitedoscillator_wrap.cxx  -lstdc++ -lm -lluajit
