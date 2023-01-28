    swig -lua -c++ -Iinclude FxDSPPanLaw.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPPanLaw.so FxDSPPanLaw_wrap.cxx -lstdc++ -lm -lluajit    
    