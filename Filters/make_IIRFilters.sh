    swig -lua -c++ -I../include -I../include IIRFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIRFilters.so IIRFilters_wrap.cxx -lstdc++ -lm -lluajit 
    