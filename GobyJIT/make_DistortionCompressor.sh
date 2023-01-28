    swig -lua -c++ -Iinclude DistortionCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DistortionCompressor.so DistortionCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    