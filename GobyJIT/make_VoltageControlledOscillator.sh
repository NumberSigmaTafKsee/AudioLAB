    swig -lua -c++ -Iinclude VoltageControlledOscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VoltageControlledOscillator.so VoltageControlledOscillator_wrap.cxx -lstdc++ -lm -lluajit    
    