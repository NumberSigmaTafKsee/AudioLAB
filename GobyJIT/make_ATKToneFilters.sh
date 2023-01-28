    swig -lua -c++ -Iinclude ATKToneFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKToneFilters.so ATKToneFilters_wrap.cxx -lstdc++ -lm -lluajit    
    