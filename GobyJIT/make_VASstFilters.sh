    swig -lua -c++ -Iinclude VASstFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASstFilters.so VASstFilters_wrap.cxx -lstdc++ -lm -lluajit    
    