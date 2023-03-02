    swig -lua -c++ -Iinclude ATKBesselFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKBesselFilters.so ATKBesselFilters_wrap.cxx -lstdc++ -lm -lluajit    
    