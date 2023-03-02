    swig -lua -c++ -Iinclude ATKRBJFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKRBJFilters.so ATKRBJFilters_wrap.cxx -lstdc++ -lm -lluajit    
    