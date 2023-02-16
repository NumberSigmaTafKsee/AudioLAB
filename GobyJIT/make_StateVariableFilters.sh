    swig -lua -c++ -Iinclude StateVariableFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o StateVariableFilters.so StateVariableFilters_wrap.cxx -lstdc++ -lm -lluajit    
    