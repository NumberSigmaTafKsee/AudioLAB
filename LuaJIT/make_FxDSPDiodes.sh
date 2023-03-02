    swig -lua -c++ -Iinclude FxDSPDiodes.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPDiodes.so FxDSPDiodes_wrap.cxx -lstdc++ -lm -lluajit    
    