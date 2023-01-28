    swig -lua -c++ -Iinclude PolyPhaseFilterBank.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PolyPhaseFilterBank.so PolyPhaseFilterBank_wrap.cxx -lstdc++ -lm -lluajit    
    