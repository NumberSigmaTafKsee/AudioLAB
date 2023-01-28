    swig -lua -c++ -Iinclude ATKUniversalFixedDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKUniversalFixedDelayLine.so ATKUniversalFixedDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    