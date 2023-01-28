    swig -lua -c++ -Iinclude FxDSPPolySaturator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPPolySaturator.so FxDSPPolySaturator_wrap.cxx -lstdc++ -lm -lluajit    
    