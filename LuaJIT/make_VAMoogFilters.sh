    swig -lua -c++ -Iinclude VAMoogFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilters.so VAMoogFilters_wrap.cxx -lstdc++ -lm -lluajit    
    