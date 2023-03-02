    swig -lua -c++ -Iinclude FxDSPMultibandFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPMultibandFilter.so FxDSPMultibandFilter_wrap.cxx -lstdc++ -lm -lluajit    
    