    swig -lua -c++ -Iinclude WDFCompressor.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WDFCompressor.h.so WDFCompressor.h_wrap.cxx -lstdc++ -lm -lluajit    
    