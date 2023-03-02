    swig -lua -c++ -Iinclude VAPolyBLEPOscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAPolyBLEPOscillator.so VAPolyBLEPOscillator_wrap.cxx -lstdc++ -lm -lluajit    
    