    swig -lua -c++ -Iinclude NoiseGenerators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o NoiseGenerators.so NoiseGenerators_wrap.cxx -lstdc++ -lm -lluajit    
    