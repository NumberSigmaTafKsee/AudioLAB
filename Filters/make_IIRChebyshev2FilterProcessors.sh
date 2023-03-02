    swig -lua -c++ -I../include -I../include IIRChebyshev2FilterProcessors.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIRChebyshev2FilterProcessors.so IIRChebyshev2FilterProcessors_wrap.cxx -lstdc++ -lm -lluajit 
    