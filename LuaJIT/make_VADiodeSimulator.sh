    swig -lua -c++ -Iinclude VADiodeSimulator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiodeSimulator.so VADiodeSimulator_wrap.cxx -lstdc++ -lm -lluajit    
    