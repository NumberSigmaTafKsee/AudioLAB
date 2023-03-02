    swig -lua -c++ -Iinclude AnalogCircuits.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AnalogCircuits.so AnalogCircuits_wrap.cxx -lstdc++ -lm -lluajit    
    