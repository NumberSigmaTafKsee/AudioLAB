    swig -lua -c++ -Iinclude FXChorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXChorus.so FXChorus_wrap.cxx -lstdc++ -lm -lluajit    
    