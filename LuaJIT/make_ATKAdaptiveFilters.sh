    swig -lua -c++ -Iinclude ATKAdaptiveFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKAdaptiveFilters.so ATKAdaptiveFilters_wrap.cxx -lstdc++ -lm -lluajit    
    