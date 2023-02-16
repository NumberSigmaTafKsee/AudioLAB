    swig -lua -c++ -Iinclude VanDerPolOscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VanDerPolOscillator.so VanDerPolOscillator_wrap.cxx -lstdc++ -lm -lluajit    
    