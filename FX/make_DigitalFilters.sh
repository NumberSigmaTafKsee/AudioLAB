    swig -lua -c++ -I../include -I../include DigitalFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o DigitalFilters.so DigitalFilters_wrap.cxx -lstdc++ -lm -lluajit 
    