    swig -lua -c++ -I../include -I../include fir_filter.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o fir_filter.h.so fir_filter.h_wrap.cxx -lstdc++ -lm -lluajit 
    