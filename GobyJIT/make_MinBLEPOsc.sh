    swig -lua -c++ -Iinclude MinBLEPOsc.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MinBLEPOsc.so MinBLEPOsc_wrap.cxx -lstdc++ -lm -lluajit    
    