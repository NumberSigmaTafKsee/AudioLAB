    swig -lua -c++ -Iinclude ATKUniversalVariableDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKUniversalVariableDelayLine.so ATKUniversalVariableDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    