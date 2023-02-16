    swig -lua -c++ -Iinclude MoogFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MoogFilters.so MoogFilters_wrap.cxx -lstdc++ -lm -lluajit    
    