    swig -lua -c++ -Iinclude DigitalFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DigitalFilters.so DigitalFilters_wrap.cxx -lstdc++ -lm -lluajit    
    