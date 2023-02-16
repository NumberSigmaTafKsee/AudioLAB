    swig -lua -c++ -Iinclude CombFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CombFilters.so CombFilters_wrap.cxx -lstdc++ -lm -lluajit    
    