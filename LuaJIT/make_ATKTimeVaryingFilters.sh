    swig -lua -c++ -Iinclude ATKTimeVaryingFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKTimeVaryingFilters.so ATKTimeVaryingFilters_wrap.cxx -lstdc++ -lm -lluajit    
    