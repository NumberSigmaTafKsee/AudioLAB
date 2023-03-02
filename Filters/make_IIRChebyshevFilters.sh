    swig -lua -c++ -I../include -I../include IIRChebyshevFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIRChebyshevFilters.so IIRChebyshevFilters_wrap.cxx -lstdc++ -lm -lluajit 
    