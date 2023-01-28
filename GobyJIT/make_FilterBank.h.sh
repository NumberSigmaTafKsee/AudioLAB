    swig -lua -c++ -Iinclude FilterBank.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FilterBank.h.so FilterBank.h_wrap.cxx -lstdc++ -lm -lluajit    
    