    swig -lua -c++ -Iinclude VAStateVariableFilter1.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStateVariableFilter1.so VAStateVariableFilter1_wrap.cxx -lstdc++ -lm -lluajit    
    