    swig -lua -c++ -Iinclude DelayFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelayFilters.so DelayFilters_wrap.cxx -lstdc++ -lm -lluajit    
    