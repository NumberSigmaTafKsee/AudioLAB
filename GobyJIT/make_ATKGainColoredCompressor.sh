    swig -lua -c++ -Iinclude ATKGainColoredCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainColoredCompressor.so ATKGainColoredCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    