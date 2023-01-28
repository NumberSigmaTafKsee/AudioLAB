    swig -lua -c++ -Iinclude WDFCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WDFCompressor.so WDFCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    