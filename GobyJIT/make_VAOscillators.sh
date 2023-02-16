    swig -lua -c++ -Iinclude VAOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAOscillators.so VAOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    