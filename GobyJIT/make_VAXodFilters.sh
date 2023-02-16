    swig -lua -c++ -Iinclude VAXodFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAXodFilters.so VAXodFilters_wrap.cxx -lstdc++ -lm -lluajit    
    