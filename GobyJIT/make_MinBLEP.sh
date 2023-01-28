    swig -lua -c++ -Iinclude MinBLEP.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MinBLEP.so MinBLEP_wrap.cxx -lstdc++ -lm -lluajit    
    