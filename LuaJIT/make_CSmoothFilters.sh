    swig -lua -c++ -Iinclude CSmoothFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CSmoothFilters.so CSmoothFilters_wrap.cxx -lstdc++ -lm -lluajit    
    