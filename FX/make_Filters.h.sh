    swig -lua -c++ -I../include -I../include Filters.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Filters.h.so Filters.h_wrap.cxx -lstdc++ -lm -lluajit 
    