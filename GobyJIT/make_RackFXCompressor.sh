    swig -lua -c++ -Iinclude RackFXCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXCompressor.so RackFXCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    