    swig -lua -c++ -Iinclude Fir1.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Fir1.h.so Fir1.h_wrap.cxx -lstdc++ -lm -lluajit    
    