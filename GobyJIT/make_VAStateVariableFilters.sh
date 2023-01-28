    swig -lua -c++ -Iinclude VAStateVariableFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStateVariableFilters.so VAStateVariableFilters_wrap.cxx -lstdc++ -lm -lluajit    
    