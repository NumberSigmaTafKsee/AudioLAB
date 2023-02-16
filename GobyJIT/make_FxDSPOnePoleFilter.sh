    swig -lua -c++ -Iinclude FxDSPOnePoleFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPOnePoleFilter.so FxDSPOnePoleFilter_wrap.cxx -lstdc++ -lm -lluajit    
    