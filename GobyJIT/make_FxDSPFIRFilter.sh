    swig -lua -c++ -Iinclude FxDSPFIRFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPFIRFilter.so FxDSPFIRFilter_wrap.cxx -lstdc++ -lm -lluajit    
    