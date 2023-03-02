    swig -lua -c++ -Iinclude BandLimitedOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BandLimitedOscillators.so BandLimitedOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    