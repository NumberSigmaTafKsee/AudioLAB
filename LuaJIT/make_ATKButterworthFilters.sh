    swig -lua -c++ -Iinclude ATKButterworthFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKButterworthFilters.so ATKButterworthFilters_wrap.cxx -lstdc++ -lm -lluajit    
    