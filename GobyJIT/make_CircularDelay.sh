    swig -lua -c++ -Iinclude CircularDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CircularDelay.so CircularDelay_wrap.cxx -lstdc++ -lm -lluajit    
    