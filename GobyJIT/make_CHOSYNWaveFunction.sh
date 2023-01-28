    swig -lua -c++ -Iinclude CHOSYNWaveFunction.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CHOSYNWaveFunction.so CHOSYNWaveFunction_wrap.cxx -lstdc++ -lm -lluajit    
    