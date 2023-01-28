    swig -lua -c++ -Iinclude FXChorus2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXChorus2.so FXChorus2_wrap.cxx -lstdc++ -lm -lluajit    
    