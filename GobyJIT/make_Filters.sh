    swig -lua -c++ -Iinclude Filters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Filters.so Filters_wrap.cxx -lstdc++ -lm -lluajit    
    