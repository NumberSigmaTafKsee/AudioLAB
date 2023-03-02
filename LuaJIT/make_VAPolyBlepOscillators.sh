    swig -lua -c++ -Iinclude VAPolyBlepOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAPolyBlepOscillators.so VAPolyBlepOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    