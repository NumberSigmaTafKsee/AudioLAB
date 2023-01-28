    swig -lua -c++ -Iinclude VAMinBlepOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMinBlepOscillators.so VAMinBlepOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    