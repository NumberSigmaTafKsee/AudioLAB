    swig -lua -c++ -Iinclude ATKTimeVaryingSVFFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKTimeVaryingSVFFilters.so ATKTimeVaryingSVFFilters_wrap.cxx -lstdc++ -lm -lluajit    
    