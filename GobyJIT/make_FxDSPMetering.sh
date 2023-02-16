    swig -lua -c++ -Iinclude FxDSPMetering.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPMetering.so FxDSPMetering_wrap.cxx -lstdc++ -lm -lluajit    
    