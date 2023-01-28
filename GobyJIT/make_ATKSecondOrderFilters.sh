    swig -lua -c++ -Iinclude ATKSecondOrderFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKSecondOrderFilters.so ATKSecondOrderFilters_wrap.cxx -lstdc++ -lm -lluajit    
    