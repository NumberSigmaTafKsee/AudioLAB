    swig -lua -c++ -Iinclude DiodeSimulator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DiodeSimulator.so DiodeSimulator_wrap.cxx -lstdc++ -lm -lluajit    
    