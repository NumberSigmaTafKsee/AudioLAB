    swig -lua -c++ -Iinclude VAStateVariableFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStateVariableFilter2.so VAStateVariableFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    