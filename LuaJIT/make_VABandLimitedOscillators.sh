    swig -lua -c++ -Iinclude VABandLimitedOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VABandLimitedOscillators.so VABandLimitedOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    