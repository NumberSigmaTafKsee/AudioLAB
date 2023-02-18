    swig -lua -c++ -I../include -I../include ResamplerProcessors.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ResamplerProcessors.so ResamplerProcessors_wrap.cxx -lstdc++ -lm -lluajit 
    