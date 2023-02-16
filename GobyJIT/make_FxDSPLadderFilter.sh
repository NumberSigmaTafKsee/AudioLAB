    swig -lua -c++ -Iinclude FxDSPLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPLadderFilter.so FxDSPLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    