    swig -lua -c++ -Iinclude PolyPhaseFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PolyPhaseFilter.so PolyPhaseFilter_wrap.cxx -lstdc++ -lm -lluajit    
    