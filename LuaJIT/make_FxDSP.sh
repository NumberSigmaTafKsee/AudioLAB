    swig -lua -c++ -Iinclude FxDSP.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSP.so FxDSP_wrap.cxx -lstdc++ -lm -lluajit    
    