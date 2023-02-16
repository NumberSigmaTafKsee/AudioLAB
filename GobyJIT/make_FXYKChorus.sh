    swig -lua -c++ -Iinclude FXYKChorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXYKChorus.so FXYKChorus_wrap.cxx -lstdc++ -lm -lluajit    
    