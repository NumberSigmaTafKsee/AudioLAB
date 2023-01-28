    swig -lua -c++ -Iinclude Noise.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Noise.h.so Noise.h_wrap.cxx -lstdc++ -lm -lluajit    
    