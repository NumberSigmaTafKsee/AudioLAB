    swig -lua -c++ -I../include -I../include IIRZolzerFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIRZolzerFilter.so IIRZolzerFilter_wrap.cxx -lstdc++ -lm -lluajit 
    