    swig -lua -c++ -Iinclude FXChorusDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXChorusDelay.so FXChorusDelay_wrap.cxx -lstdc++ -lm -lluajit    
    