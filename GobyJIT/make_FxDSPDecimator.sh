    swig -lua -c++ -Iinclude FxDSPDecimator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPDecimator.so FxDSPDecimator_wrap.cxx -lstdc++ -lm -lluajit    
    