    swig -lua -c++ -Iinclude FxDSPFFT.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPFFT.so FxDSPFFT_wrap.cxx -lstdc++ -lm -lluajit    
    