    swig -lua -c++ -Iinclude VAAnalogFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAAnalogFilters.so VAAnalogFilters_wrap.cxx -lstdc++ -lm -lluajit    
    