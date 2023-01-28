    swig -lua -c++ -Iinclude AllPassFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AllPassFilters.so AllPassFilters_wrap.cxx -lstdc++ -lm -lluajit    
    