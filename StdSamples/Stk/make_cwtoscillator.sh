swig -lua -c++ -Iinclude cwtoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cwtoscillator.so cwtoscillator_wrap.cxx  -lstdc++ -lm -lluajit
