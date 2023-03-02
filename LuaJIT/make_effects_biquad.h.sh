    swig -lua -c++ -Iinclude effects_biquad.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o effects_biquad.h.so effects_biquad.h_wrap.cxx -lstdc++ -lm -lluajit    
    