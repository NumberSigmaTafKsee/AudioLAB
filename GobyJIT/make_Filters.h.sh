    swig -lua -c++ -Iinclude Filters.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Filters.h.so Filters.h_wrap.cxx -lstdc++ -lm -lluajit    
    