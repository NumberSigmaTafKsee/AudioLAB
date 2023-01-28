    swig -lua -c++ -Iinclude FxDSPBiquads.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPBiquads.so FxDSPBiquads_wrap.cxx -lstdc++ -lm -lluajit    
    