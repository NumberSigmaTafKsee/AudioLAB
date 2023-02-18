    swig -lua -c++ -I../include -I../include Noise.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Noise.h.so Noise.h_wrap.cxx -lstdc++ -lm -lluajit 
    