    swig -lua -c++ -Iinclude VAStateVariableCombFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStateVariableCombFilter.so VAStateVariableCombFilter_wrap.cxx -lstdc++ -lm -lluajit    
    