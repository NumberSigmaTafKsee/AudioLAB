    swig -lua -c++ -Iinclude CTAGDRCCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CTAGDRCCompressor.so CTAGDRCCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    