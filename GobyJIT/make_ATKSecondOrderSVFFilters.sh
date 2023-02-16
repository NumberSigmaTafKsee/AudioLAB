    swig -lua -c++ -Iinclude ATKSecondOrderSVFFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKSecondOrderSVFFilters.so ATKSecondOrderSVFFilters_wrap.cxx -lstdc++ -lm -lluajit    
    