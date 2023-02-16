    swig -lua -c++ -Iinclude VAWDFCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAWDFCompressor.so VAWDFCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    