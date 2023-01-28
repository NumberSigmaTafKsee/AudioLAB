    swig -lua -c++ -Iinclude Waveshaping.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Waveshaping.so Waveshaping_wrap.cxx -lstdc++ -lm -lluajit    
    