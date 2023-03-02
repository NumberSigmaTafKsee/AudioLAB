    swig -lua -c++ -Iinclude effects_biquad.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o effects_biquad.cpp.so effects_biquad.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    