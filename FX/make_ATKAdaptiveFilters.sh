    swig -lua -c++ -I../include -I../include ATKAdaptiveFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKAdaptiveFilters.so ATKAdaptiveFilters_wrap.cxx -lstdc++ -lm -lluajit 
    