    swig -lua -c++ -Iinclude SstFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SstFilters.so SstFilters_wrap.cxx -lstdc++ -lm -lluajit    
    