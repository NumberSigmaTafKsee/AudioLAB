    swig -lua -c++ -Iinclude ATKGainCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainCompressor.so ATKGainCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    