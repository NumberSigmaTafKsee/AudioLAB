    swig -lua -c++ -Iinclude Waveshapers.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Waveshapers.so Waveshapers_wrap.cxx -lstdc++ -lm -lluajit    
    