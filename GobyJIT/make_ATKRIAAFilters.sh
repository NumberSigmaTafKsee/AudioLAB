    swig -lua -c++ -Iinclude ATKRIAAFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKRIAAFilters.so ATKRIAAFilters_wrap.cxx -lstdc++ -lm -lluajit    
    