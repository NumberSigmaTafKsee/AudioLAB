    swig -lua -c++ -Iinclude FxDSPProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPProcessors.so FxDSPProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    