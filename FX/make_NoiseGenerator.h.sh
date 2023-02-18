    swig -lua -c++ -I../include -I../include NoiseGenerator.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o NoiseGenerator.h.so NoiseGenerator.h_wrap.cxx -lstdc++ -lm -lluajit 
    