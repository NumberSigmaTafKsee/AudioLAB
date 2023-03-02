    swig -lua -c++ -I../include -I../include IIREllipticFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIREllipticFilter.so IIREllipticFilter_wrap.cxx -lstdc++ -lm -lluajit 
    