    swig -lua -c++ -Iinclude NoiseGenerator.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o NoiseGenerator.h.so NoiseGenerator.h_wrap.cxx -lstdc++ -lm -lluajit    
    