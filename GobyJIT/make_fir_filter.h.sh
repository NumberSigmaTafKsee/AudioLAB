    swig -lua -c++ -Iinclude fir_filter.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fir_filter.h.so fir_filter.h_wrap.cxx -lstdc++ -lm -lluajit    
    