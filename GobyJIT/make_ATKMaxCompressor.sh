    swig -lua -c++ -Iinclude ATKMaxCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKMaxCompressor.so ATKMaxCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    