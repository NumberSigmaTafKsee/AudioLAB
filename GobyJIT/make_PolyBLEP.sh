    swig -lua -c++ -Iinclude PolyBLEP.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PolyBLEP.so PolyBLEP_wrap.cxx -lstdc++ -lm -lluajit    
    