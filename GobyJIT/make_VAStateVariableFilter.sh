    swig -lua -c++ -Iinclude VAStateVariableFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStateVariableFilter.so VAStateVariableFilter_wrap.cxx -lstdc++ -lm -lluajit    
    